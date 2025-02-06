{-# LANGUAGE OverloadedStrings #-}

module Main
  ( main
  ) where

import           Data.ByteString.Lazy.Char8 (pack)
import           Data.Either                (isLeft)
import qualified Data.Map                   as Map
import           Data.Set                   (singleton)
import qualified Data.Set                   as Set
import           Exploration
import           Metabolome
import           Test.Hspec

spec :: Spec
spec = do
  let r1 = ("R1", Reaction (Set.fromList ["A"]) (Set.fromList ["B"]) "A -> B")
      r1Rev =
        ("R1_r", Reaction (Set.fromList ["B"]) (Set.fromList ["A"]) "A -> B")
      r2 = ("R2", Reaction (Set.fromList ["B"]) (Set.fromList ["C"]) "B -> C")
      r3 = ("R3", Reaction (Set.fromList ["B"]) (Set.fromList ["C"]) "B -> C")
      r3' = ("R3", Reaction (Set.fromList ["X"]) (Set.fromList ["D"]) "X -> D")
      r4 = ("R4", Reaction (Set.fromList ["A", "X"]) (Set.fromList ["B"]) "")
  describe "readReactions" $ do
    it "parses an empty JSON" $ do
      readReactions (pack "{}") `shouldBe` Right Map.empty
    it "parses a valid JSON into a ReactionMap" $ do
      let jsonInput =
            "{ \"R1\":{\"Reactants\":{\"A\":[\"metaboliteA\"]},\"Products\":{\"B\":[\"metaboliteB\"]},\"Equation\":\"A -> B\"}}"
      readReactions (pack jsonInput) `shouldBe` Right (Map.fromList [r1])
    it "fails on invalid JSON" $ do
      let jsonInput = "[ \"R\" ]"
          result = readReactions (pack jsonInput)
      result `shouldSatisfy` isLeft
  describe "mirrorReactions" $ do
    it "creates reversed reactions with a modified key" $ do
      let reactions = Map.fromList [r1]
          mirrored = mirrorReactions reactions
          expected = Map.fromList [r1Rev]
      mirrored `shouldBe` expected
  describe "Exploration functions" $ do
    it "computes the correct tree size" $ do
      let tree =
            Tree
              Nothing
              (Set.fromList ["A"])
              (Set.fromList ["A", "B"])
              (Just [])
              Set.empty
      treeSize tree `shouldBe` 1
    it "computes all products in a tree" $ do
      let tree =
            Tree
              Nothing
              (Set.fromList ["A"])
              (Set.fromList ["A", "B"])
              (Just
                 [ Tree
                     (Just "R1")
                     (Set.fromList ["C"])
                     (Set.fromList ["A", "B", "C"])
                     Nothing
                     Set.empty
                 ])
              Set.empty
      allProducts tree `shouldBe` Set.fromList ["A", "C"]
    it "validates strict policy correctly" $ do
      let tree =
            Tree
              Nothing
              (Set.fromList ["A"])
              (Set.fromList ["A"])
              Nothing
              Set.empty
      strict tree (snd r1) `shouldBe` True
    it "validates permissive policy correctly" $ do
      let reaction = Reaction (Set.fromList ["A"]) (Set.fromList ["C"]) "A -> C"
          tree =
            Tree
              Nothing
              (Set.fromList ["A"])
              (Set.fromList ["A"])
              Nothing
              Set.empty
      permissive tree reaction `shouldBe` True
    it "expands network correctly with depth" $ do
      let reactions = Map.fromList [r1, r2]
          metabolites = Set.fromList ["A"]
          expandedTree =
            expansion metabolites reactions (singleton "A") strict 2
      treeSize expandedTree `shouldBe` 3
    it "expands network correctly with redundant reaction" $ do
      let reactions = Map.fromList [r1, r2, r3]
          metabolites = Set.fromList ["A"]
          expandedTree =
            expansion metabolites reactions (singleton "A") strict 2
      treeSize expandedTree `shouldBe` 4
    it "expands network correctly disconnected reaction" $ do
      let reactions = Map.fromList [r1, r2, r3']
          metabolites = Set.fromList ["A"]
          expandedTree =
            expansion metabolites reactions (singleton "A") strict 2
      treeSize expandedTree `shouldBe` 3
    it "expands network correctly permissively" $ do
      let reactions = Map.fromList [r4]
          metabolites = Set.fromList ["A"]
          expandedTree =
            expansion metabolites reactions (singleton "A") permissive 2
      treeSize expandedTree `shouldBe` 2
    it "not expands known cmpds permissively" $ do
      let reactions = Map.fromList [r4]
          metabolites = Set.fromList ["A", "B"]
          expandedTree =
            expansion metabolites reactions (singleton "A") permissive 2
      treeSize expandedTree `shouldBe` 1
    it "not expands network correctly permissively" $ do
      let reactions = Map.fromList [r2]
          metabolites = Set.fromList ["A"]
          expandedTree =
            expansion metabolites reactions (singleton "A") permissive 2
      treeSize expandedTree `shouldBe` 1

main :: IO ()
main = hspec spec
