{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE LambdaCase            #-}
{-# LANGUAGE OverloadedStrings     #-}
{-# LANGUAGE ScopedTypeVariables   #-}

module Main
  ( main
  ) where

import           Control.Monad           (forM_, when)
import           Data.Aeson              (eitherDecode, encodeFile)
import qualified Data.ByteString.Lazy    as BS
import qualified Data.Configurator       as Cfg
import           Data.Configurator.Types (Config)
import           Data.Csv
import           Data.Foldable           (toList)
import qualified Data.Map                as Map
import           Data.Set                ((\\))
import qualified Data.Set                as Set
import           Data.Text               (Text)
import           System.Environment      (getArgs)
import           System.Exit             (ExitCode (..), exitWith)

import           Draw
import           Exploration
import           Metabolome

readPolicy :: Config -> IO Policy
readPolicy config =
  Cfg.require config "exploration.policy" >>= \case
    "strict" -> pure strict
    "permissive" -> pure permissive
    "relevant" -> pure relevant
    (p :: String) -> do
      putStrLn
        $ "Error: Unknown policy "
            <> show p
            <> ". Use \"strict\" or \"permissive\"."
      exitWith (ExitFailure 1)

readPresentMetabolites :: Config -> IO MetaboliteNames
readPresentMetabolites config =
  Cfg.lookup config "input.model" >>= \case
    Nothing -> return Set.empty
    Just filepath -> do
      contentJSON <- BS.readFile filepath
      case eitherDecode contentJSON of
        Left err -> do
          putStrLn $ "Warning: Error decoding input.model: " <> err
          return Set.empty
        Right content -> return (modelMetabolites content)

readReactome :: Config -> IO ReactionMap
readReactome config = do
  reactomeJSON <- Cfg.require config "input.reactome" >>= BS.readFile
  case eitherDecode reactomeJSON of
    Left err -> do
      putStrLn $ "Error: " <> err
      exitWith (ExitFailure 1)
    Right reactome -> return $ reactome <> mirrorReactions reactome

expandTree :: Config -> IO ()
expandTree config = do
  putStrLn "-> Reading data"
  reactome <- readReactome config
  metabolites <- readPresentMetabolites config
  policy <- readPolicy config
  depth <- Cfg.require config "exploration.depth"
  initial <- Set.fromList <$> Cfg.require config "exploration.initial"
  putStrLn
    $ "-> BFS to depth " <> show depth <> "; starting at " <> show initial
  let expandedTree = expansion metabolites reactome initial policy depth
  let newCompounds = allProducts expandedTree \\ metabolites
  putStrLn $ "   Search tree size: " <> show (treeSize expandedTree)
  putStrLn $ "   New compounds: " <> show (length newCompounds)
  writeOutput config reactome expandedTree

data CompoundRecord = CompoundRecord
  { keggId   :: MetaboliteId
  , cmpdName :: Text
  } deriving (Show)

instance FromNamedRecord CompoundRecord where
  parseNamedRecord r = CompoundRecord <$> r .: "KEGG" <*> r .: "name"

readCompoundMap :: FilePath -> IO (Map.Map MetaboliteId Text)
readCompoundMap filePath = do
  decodeByName <$> BS.readFile filePath >>= \case
    Left err -> do
      putStrLn $ "Error parsing CSV: " <> err
      return Map.empty
    Right (_, records) -> do
      putStrLn $ "   Compound names read from " <> filePath
      return $ Map.fromList [(keggId r, cmpdName r) | r <- toList records]

data CompoundSimilarity = CompoundSimilarity
  { kegg       :: MetaboliteId
  , similarity :: Float
  } deriving (Show)

instance FromNamedRecord CompoundSimilarity where
  parseNamedRecord r =
    CompoundSimilarity <$> r .: "keggId" <*> r .: "similarity"

readCompoundSimilarities :: FilePath -> IO (Map.Map MetaboliteId Float)
readCompoundSimilarities filePath = do
  decodeByName <$> BS.readFile filePath >>= \case
    Left err -> do
      putStrLn $ "Error parsing CSV: " <> err
      return Map.empty
    Right (_, records) -> do
      putStrLn $ "   Compound similarities read form " <> filePath
      return $ Map.fromList [(kegg r, similarity r) | r <- toList records]

writeGraph :: Config -> ReactionMap -> Tree -> FilePath -> IO ()
writeGraph config reactions tree file = do
  names <-
    Cfg.lookup config "input.compound_info"
      >>= maybe (pure Map.empty) readCompoundMap
  similarities <-
    Cfg.lookup config "input.similarities"
      >>= maybe (pure Map.empty) readCompoundSimilarities
  Cfg.lookup config "exploration.ascii"
    >>= mapM_ (flip when (putStr (asciiTree names tree)))
  BS.writeFile file $ plotTreeAsGraph reactions names similarities tree
  putStrLn $ "   Graph written to " <> file

writeExplorationJSON :: Tree -> FilePath -> IO ()
writeExplorationJSON tree file = do
  encodeFile file tree
  putStrLn $ "   Exploration data written to " <> file

writeOutput :: Config -> ReactionMap -> Tree -> IO ()
writeOutput config reactome tree = do
  Cfg.lookup config "exploration.output" >>= mapM_ (writeExplorationJSON tree)
  Cfg.lookup config "exploration.graph"
    >>= mapM_ (writeGraph config reactome tree)

main :: IO ()
main = do
  args <- getArgs
  when (null args) $ do
    putStrLn "Error: Supply a configuration file."
    exitWith (ExitFailure 1)
  forM_ args $ \configFile -> do
    putStrLn $ "-> Configuration: " <> configFile
    expandTree =<< Cfg.load [Cfg.Required configFile]
