{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Main
  ( main
  ) where

import           Data.Aeson
import qualified Data.ByteString.Lazy as BS
import qualified Data.Configurator    as Cfg
import           Data.Set             ((\\))

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

expandTree :: IO ()
expandTree = do
  cfgFile <- Cfg.load ["config.cfg"]
  putStrLn "Reading data"
  reactomeJSON <- BS.readFile =<< Cfg.require cfgFile "input.reactome"
  metabolismJSON <- BS.readFile =<< Cfg.require cfgFile "input.model"
  (policyPar :: String) <- Cfg.require cfgFile "exploration.policy"
  let policy =
        case policyPar of
          "strict" -> strict
          "permissive" -> permissive
          _ ->
            error
              $ "unknown policy " <> policyPar <> ". Use strict or permissive."
  case readData reactomeJSON metabolismJSON of
    Left err -> putStrLn $ "Error: " ++ err
    Right (metabolites, reactome) -> do
      depth <- Cfg.require cfgFile "exploration.depth"
      initial <- Cfg.require cfgFile "exploration.initial"
      putStrLn
        $ "BFS to depth " <> show depth <> "; starting at " <> show initial
      let expandedTree = expansion metabolites reactome initial policy depth
      let newCompounds = allProducts expandedTree \\ metabolites
      putStrLn $ "Search tree size: " <> show (treeSize expandedTree)
      putStrLn $ "New compounds: " <> show (length newCompounds)
      searchSpaceFilePath <- Cfg.require cfgFile "exploration.output"
      putStrLn $ "Search tree written to " <> searchSpaceFilePath
      encodeFile searchSpaceFilePath expandedTree

main :: IO ()
main = expandTree
