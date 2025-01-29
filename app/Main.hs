{-# LANGUAGE LambdaCase          #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Main
  ( main
  ) where

import           Control.Monad           (forM_, when)
import           Data.Aeson              (encodeFile)
import qualified Data.ByteString.Lazy    as BS
import qualified Data.Configurator       as Cfg
import           Data.Configurator.Types (Config)
import           Data.Csv
import           Data.Foldable           (toList)
import qualified Data.Map                as Map
import           Data.Set                ((\\))
import           Data.Text               (Text)
import           System.Environment      (getArgs)
import           System.Exit             (ExitCode (..), exitWith)

import           Draw
import           Exploration
import           Metabolome

readData ::
     BS.ByteString
  -> BS.ByteString
  -> Either String (MetaboliteNames, ReactionMap)
readData reactomeJSON metabolismJSON = do
  (_, reactome) <- readReactions reactomeJSON
  let reactomeRev = mirrorReactions reactome
  (modelMetabolites, _) <- readReactions metabolismJSON
  return (modelMetabolites, reactome <> reactomeRev)

expandTree :: Config -> IO ()
expandTree config = do
  putStrLn "-> Reading data"
  reactomeJSON <- BS.readFile =<< Cfg.require config "input.reactome"
  metabolismJSON <- BS.readFile =<< Cfg.require config "input.model"
  policy <-
    Cfg.require config "exploration.policy" >>= \case
      "strict" -> pure strict
      "permissive" -> pure permissive
      (p :: String) -> do
        putStrLn
          $ "Error: Unknown policy "
              <> show p
              <> ". Use \"strict\" or \"permissive\"."
        exitWith (ExitFailure 1)
  case readData reactomeJSON metabolismJSON of
    Left err -> putStrLn $ "Error: " ++ err
    Right (metabolites, reactome) -> do
      depth <- Cfg.require config "exploration.depth"
      initial <- Cfg.require config "exploration.initial"
      putStrLn
        $ "-> BFS to depth " <> show depth <> "; starting at " <> show initial
      let expandedTree = expansion metabolites reactome initial policy depth
      let newCompounds = allProducts expandedTree \\ metabolites
      putStrLn $ "   Search tree size: " <> show (treeSize expandedTree)
      putStrLn $ "   New compounds: " <> show (length newCompounds)
      writeOutput config expandedTree

data CompoundRecord = CompoundRecord
  { cmpdName :: Text
  , keggId   :: MetaboliteId
  } deriving (Show)

instance FromNamedRecord CompoundRecord where
  parseNamedRecord r = CompoundRecord <$> r .: "name" <*> r .: "KEGG"

readCompoundMap :: FilePath -> IO (Map.Map MetaboliteId Text)
readCompoundMap filePath = do
  decodeByName <$> BS.readFile filePath >>= \case
    Left err -> do
      putStrLn $ "Error parsing CSV: " ++ err
      return Map.empty
    Right (_, records) -> do
      putStrLn $ "   Compound names read from " <> filePath
      return $ Map.fromList [(keggId r, cmpdName r) | r <- toList records]

writeGraph :: Config -> Tree -> FilePath -> IO ()
writeGraph config tree file = do
  names <-
    Cfg.lookup config "input.compound_info"
      >>= maybe (pure Map.empty) readCompoundMap
  BS.writeFile file $ plotTreeAsGraph names tree
  putStrLn $ "   Graph written to " <> file

writeExplorationJSON :: Tree -> FilePath -> IO ()
writeExplorationJSON tree file = do
  encodeFile file tree
  putStrLn $ "   Exploration data written to " <> file

writeOutput :: Config -> Tree -> IO ()
writeOutput config tree = do
  Cfg.lookup config "exploration.output" >>= mapM_ (writeExplorationJSON tree)
  Cfg.lookup config "exploration.graph" >>= mapM_ (writeGraph config tree)

main :: IO ()
main = do
  args <- getArgs
  when (null args) $ do
    putStrLn "Error: Supply a configuration file."
    exitWith (ExitFailure 1)
  forM_ args $ \configFile -> do
    putStrLn $ "-> Configuration: " <> configFile
    expandTree =<< Cfg.load [Cfg.Required configFile]
