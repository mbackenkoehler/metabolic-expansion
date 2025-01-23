{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Exploration
  ( expansion
  , allProducts
  , treeSize
  , strict
  , permissive
  ) where

import           Data.Aeson
import           Data.Hashable (hash)
import           Data.List     (sort)
import qualified Data.Map      as Map
import           Data.Set      (Set, intersection, isSubsetOf, singleton, union,
                                (\\))
import qualified Data.Set      as Set
import           Metabolome

data Tree = Tree
  { incoming :: Maybe ReactionId
  , novel    :: Set MetaboliteId
  , present  :: Set MetaboliteId
  , children :: Maybe [Tree]
  } deriving (Show)

instance ToJSON Tree where
  toJSON tree =
    object
      [ "hash" .= nodeHash tree
      , "incoming" .= tree.incoming
      , "novel" .= tree.novel
      , "present" .= tree.present
      , "children" .= tree.children
      ]

allProducts :: Tree -> Set MetaboliteId
allProducts tree =
  Set.union tree.novel
    $ maybe Set.empty (Set.unions . fmap allProducts) (tree.children)

treeSize :: Tree -> Int
treeSize = (+ 1) . maybe 0 (sum . fmap treeSize) . children

nodeHash :: Tree -> Int
nodeHash node = hash (hashSet node.present, hashSet node.novel)
  where
    hashSet = hash . mconcat . sort . Set.toList

type Policy = Reaction -> Tree -> Bool

strict :: Policy
strict reaction node = stoichValidity reaction node && interesting reaction node

stoichValidity :: Policy
stoichValidity reaction node = reaction.reactants `isSubsetOf` node.present

interesting :: Policy
interesting reaction node =
  not (null (reaction.reactants `intersection` node.novel))
    && not (reaction.products `isSubsetOf` node.present)

permissive :: Policy
permissive = interesting

expand :: ReactionMap -> Policy -> Tree -> Tree
expand reactions policy node =
  Tree
    { incoming = node.incoming
    , novel = node.novel
    , present = node.present
    , children =
        case node.children of
          Just nodes -> Just (expand reactions policy <$> nodes)
          Nothing ->
            Just
              [ Tree
                { incoming = Just reactionId
                , novel = reaction.products \\ node.present
                , present = reaction.products `union` node.present
                , children = Nothing
                }
              | (reactionId, reaction) <- Map.assocs reactions
              , policy reaction node
              ]
    }

expansion ::
     Set MetaboliteId -> ReactionMap -> MetaboliteId -> Policy -> Int -> Tree
expansion metabolites reactions startingMetabolite policy depth =
  let startSet = singleton startingMetabolite
      root =
        Tree
          { incoming = Nothing
          , novel = startSet
          , present = metabolites `union` startSet
          , children = Nothing
          }
   in iterate (expand reactions policy) root !! depth
