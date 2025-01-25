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

readCompoundMap :: FilePath -> IO (Map.Map Text Text)
readCompoundMap filePath = do
  csvData <- BS.readFile filePath
  case decodeByName csvData of
    Left err -> do
      putStrLn $ "Error parsing CSV: " ++ err
      return Map.empty
    Right (_, records) -> do
      putStrLn $ "   Compound names read from " <> filePath
      return $ Map.fromList [(keggId rec, cmpdName rec) | rec <- toList records]

readCompoundInfo :: Config -> IO (Map.Map MetaboliteId Text)
readCompoundInfo config = do
  Cfg.lookup config "input.compound_info" >>= \case
    Nothing -> return Map.empty
    Just filename -> readCompoundMap filename

writeOutput :: Config -> Tree -> IO ()
writeOutput config tree = do
  searchSpaceFilePath <- Cfg.require config "exploration.output"
  putStrLn $ "   Search tree written to " <> searchSpaceFilePath
  encodeFile searchSpaceFilePath tree
  Cfg.lookup config "exploration.graph" >>= \case
    Nothing -> return ()
    Just outputFile -> do
      names <- readCompoundInfo config
      BS.writeFile outputFile $ plotTreeAsGraph names tree
      putStrLn $ "   Graph written to " <> outputFile

main :: IO ()
main = do
  args <- getArgs
  when (null args) $ do
    putStrLn "Error: Supply a configuration file."
    exitWith (ExitFailure 1)
  forM_ args $ \configFile -> do
    putStrLn $ "-> Configuration: " <> configFile
    expandTree =<< Cfg.load [Cfg.Required configFile]
