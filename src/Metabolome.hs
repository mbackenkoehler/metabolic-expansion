{-# LANGUAGE DeriveGeneric       #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE RecordWildCards     #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Metabolome
  ( readReactions
  , mirrorReactions
  , MetaboliteId
  , MetaboliteNames
  , ReactionMap
  , Reaction(..)
  , ReactionId
  ) where

import           Data.Aeson
import           Data.Aeson.Types
import           Data.ByteString.Lazy (ByteString)
import           Data.Map             (Map, assocs)
import qualified Data.Map             as Map
import           Data.Set             (Set)
import qualified Data.Set             as Set
import           Data.Text            (Text)
import           GHC.Generics

type ReactionId = Text

type MetaboliteId = Text

type MetaboliteNames = Set MetaboliteId

data Reaction = Reaction
  { reactants :: Set MetaboliteId
  , products  :: Set MetaboliteId
  , equation  :: Maybe Text -- contains stoichiometry
  } deriving (Show, Generic)

instance FromJSON Reaction where
  parseJSON =
    withObject "Reaction" $ \v ->
      Reaction
        <$> (Set.fromList . Map.keys
               <$> (v .: "Reactants" :: Parser (Map MetaboliteId [Text])))
        <*> (Set.fromList . Map.keys
               <$> (v .: "Products" :: Parser (Map MetaboliteId [Text])))
        <*> v .: "Equation_cid"

type ReactionMap = Map ReactionId Reaction

readReactions :: ByteString -> Either String (MetaboliteNames, ReactionMap)
readReactions jsonInput = do
  reactions <- eitherDecode jsonInput
  let metabolites = buildMetaboliteList reactions
  return (metabolites, reactions)
  where
    buildMetaboliteList reactions =
      let combine reaction = Set.union (reactants reaction) (products reaction)
       in Set.unions (combine <$> reactions)

mirrorReactions :: ReactionMap -> ReactionMap
mirrorReactions =
  Map.fromList . fmap (\(ident, r) -> (ident <> "_r", flipReaction r)) . assocs
  where
    flipReaction Reaction {reactants = rs, products = ps, ..} =
      Reaction {reactants = ps, products = rs, ..}
