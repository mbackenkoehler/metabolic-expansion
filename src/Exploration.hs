{-# LANGUAGE DeriveAnyClass      #-}
{-# LANGUAGE DeriveGeneric       #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Exploration
  ( expansion
  , allProducts
  , treeSize
  , strict
  , permissive
  , relevant
  , nodeHash
  , Policy
  , Tree(..)
  , MetaboliteId
  ) where

import           Control.Parallel.Strategies
import           Data.Aeson
import           Data.Hashable               (hash)
import           Data.List                   (sort)
import qualified Data.Map                    as Map
import           Data.Set                    (Set, intersection, isSubsetOf,
                                              union, (\\))
import qualified Data.Set                    as Set
import           GHC.Generics                (Generic)
import           Metabolome

data Tree = Tree
  { incoming :: Maybe ReactionId
  , novel    :: Set MetaboliteId
  , present  :: Set MetaboliteId
  , children :: Maybe [Tree]
  , missing  :: Set MetaboliteId
  } deriving (Show, Generic, NFData)

instance ToJSON Tree where
  toJSON tree =
    object
      [ "hash" .= nodeHash tree
      , "incoming" .= tree.incoming
      , "novel" .= tree.novel
      , "present" .= tree.present
      , "children" .= tree.children
      , "missing" .= tree.missing
      ]

allProducts :: Tree -> Set MetaboliteId
allProducts tree =
  Set.union tree.novel
    $ maybe Set.empty (Set.unions . fmap allProducts) tree.children

treeSize :: Tree -> Int
treeSize = (+ 1) . maybe 0 (sum . fmap treeSize) . children

nodeHash :: Tree -> Int
nodeHash node =
  hash [hashSet node.present, hashSet node.novel, hashSet node.missing]
  where
    hashSet = hash . mconcat . sort . Set.toList

type Policy = Tree -> Reaction -> Bool

strict :: Policy
strict n r = stoichiometry n r && relevant n r && novelty n r

stoichiometry :: Policy
stoichiometry n r = r.reactants `isSubsetOf` n.present

relevant :: Policy
relevant n r = not (null (r.reactants `intersection` n.novel))

novelty :: Policy
novelty n r = not (r.products `isSubsetOf` n.present)

permissive :: Policy
permissive n r =
  relevant n r && not (null (r.products \\ n.present \\ r.reactants))

expand :: ReactionMap -> Policy -> Tree -> Tree
expand reactions policy node =
  node
    { children =
        case node.children of
          Just nodes -> Just (parMap rdeepseq (expand reactions policy) nodes)
          Nothing    -> Just expandedChildren
    }
  where
    expandedChildren =
      newChild <$> filter (policy node . snd) (Map.assocs reactions)
    newChild (reactionId, reaction :: Reaction) =
      Tree
        { incoming = Just reactionId
        , novel = reaction.products \\ node.present
        , present = reaction.products `union` node.present
        , children = Nothing
        , missing = reaction.reactants \\ node.present
        }

expansion ::
     Set MetaboliteId
  -> ReactionMap
  -> Set MetaboliteId
  -> Policy
  -> Int
  -> Tree
expansion metabolites reactions startSet policy depth =
  let root =
        Tree
          { incoming = Nothing
          , novel = startSet
          , present = metabolites `union` startSet
          , children = Nothing
          , missing = Set.empty
          }
   in iterate (expand reactions policy) root !! depth
