{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Products
  ( expansion
  , allProducts
  , treeSize
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
      , "novel" .= novel tree
      , "children" .= children tree
      ]

allProducts :: Tree -> Set MetaboliteId
allProducts tree =
  Set.union tree.novel
    $ maybe Set.empty (Set.unions . fmap allProducts) (children tree)

treeSize :: Tree -> Int
treeSize = (+ 1) . maybe 0 (sum . fmap treeSize) . children

nodeHash :: Tree -> Int
nodeHash node = hash (hashSet node.present, hashSet node.novel)
  where
    hashSet = hash . mconcat . sort . Set.toList

expand :: ReactionMap -> Tree -> Tree
expand reactions node =
  Tree
    { incoming = node.incoming
    , novel = node.novel
    , present = node.present
    , children =
        case node.children of
          Just nodes -> Just (expand reactions <$> nodes)
          Nothing ->
            Just
              [ Tree
                { incoming = Just reactionId
                , novel = products reaction \\ present node
                , present = products reaction `union` present node
                , children = Nothing
                }
              | (reactionId, reaction) <- Map.assocs reactions
              , reactants reaction `isSubsetOf` node.present
              , not (null (reaction.reactants `intersection` node.novel))
                  && not (reaction.products `isSubsetOf` node.present)
              ]
    }

expansion :: Set MetaboliteId -> ReactionMap -> MetaboliteId -> Int -> Tree
expansion metabolites reactions startingMetabolite depth =
  let root =
        Tree
          { incoming = Nothing
          , novel = singleton startingMetabolite
          , present = metabolites
          , children = Nothing
          }
   in iterate (expand reactions) root !! depth
