{-# LANGUAGE DeriveAnyClass        #-}
{-# LANGUAGE DeriveGeneric         #-}
{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE OverloadedRecordDot   #-}
{-# LANGUAGE OverloadedStrings     #-}

module Draw
  ( plotTreeAsGraph
  ) where

import           Data.Aeson
import qualified Data.ByteString.Lazy as BL
import           Data.Hashable        (hash)
import           Data.List            (nub)
import           Data.Map             (Map, (!?))
import qualified Data.Map             as Map
import           Data.Maybe           (fromMaybe, isJust)
import qualified Data.Set             as Set
import           Data.Text
import qualified Data.Text.Encoding   as TE
import           GHC.Generics
import           Lucid
import           Lucid.Base

import           Exploration
import           Metabolome

type NodeId = Int

data Node = Node
  { color :: Text
  , size  :: Float
  , id    :: NodeId
  , label :: Text
  , shape :: Text
  } deriving (Generic, Show, Eq, ToJSON)

data Edge = Edge
  { arrows :: Text
  , title  :: Text
  , from   :: NodeId
  , to     :: NodeId
  } deriving (Generic, Show, Eq, ToJSON)

type Graph = ([Node], [Edge])

plotTreeAsGraph :: ReactionMap -> Map MetaboliteId Text -> Tree -> BL.ByteString
plotTreeAsGraph reactions names =
  renderBS . pyvis . extractGraph reactions names

extractGraph :: ReactionMap -> Map MetaboliteId Text -> Tree -> Graph
extractGraph reactions names tree = (Map.elems nodeMap, nub edges'')
  where
    (nodeMap, edges'') = extractGraph' 40 tree
    extractGraph' :: Float -> Tree -> (Map NodeId Node, [Edge])
    extractGraph' nodeSize t =
      let curHash = nodeHash t
          node =
            Node
              { color =
                  if isJust t.incoming
                    then "blue"
                    else "green"
              , size = nodeSize
              , id = curHash
              , label =
                  intercalate ", "
                    $ [fromMaybe c (names !? c) | c <- Set.toList t.novel]
              , shape = "dot"
              }
          missingNodes =
            Map.fromList
              [ ( hash m
                , Node
                    { color = "red"
                    , size = nodeSize
                    , id = hash m
                    , label = fromMaybe m (names !? m)
                    , shape = "dot"
                    })
              | m <- Set.toList t.missing
              ]
          reactionEq = maybe "" equation . (=<<) (reactions !?)
          missingEdges =
            [ Edge
              { arrows = "to"
              , from = hash m
              , to = curHash
              , title = reactionEq t.incoming
              }
            | m <- Set.toList t.missing
            ]
          descendants = fromMaybe [] t.children
          edges =
            [ Edge
              { arrows = "to"
              , from = curHash
              , to = nodeHash c
              , title = reactionEq c.incoming
              }
            | c <- descendants
            ]
          newNodeSize = nodeSize * 0.8
          (nodes', edges') = unzip $ extractGraph' newNodeSize <$> descendants
       in ( Map.insert curHash node (Map.unions (missingNodes : nodes'))
          , Prelude.concat (missingEdges : edges : edges'))

referrerpolicy_ :: Text -> Attribute
referrerpolicy_ = makeAttribute "referrerpolicy"

cssStyle :: Text
cssStyle =
  "#mynetwork {width: 100%; height: 1200px; background-color: #ffffff; border: 1px solid lightgray; position: relative; float: left;} "
    <> "#loadingBar {position: absolute; top: 0; left: 0; width: 100%; height: 1200px; background-color: rgba(200,200,200,0.8); transition: all 0.5s ease; opacity: 1;} "
    <> "#bar {position: absolute; top: 0; left: 0; width: 20px; height: 20px; margin: auto; border-radius: 11px; border: 2px solid rgba(30,30,30,0.05); background: rgb(0,173,246); box-shadow: 2px 0px 4px rgba(0,0,0,0.4);} "
    <> "#border {position: absolute; top: 10px; left: 10px; width: 500px; height: 23px; margin: auto; box-shadow: 0px 0px 4px rgba(0,0,0,0.2); border-radius: 10px;} "
    <> "#text {position: absolute; top: 8px; left: 530px; width: 30px; height: 50px; margin: auto; font-size: 22px; color: #000000;} "
    <> "div.outerBorder {position: relative; top: 400px; width: 600px; height: 44px; margin: auto; border: 8px solid rgba(0,0,0,0.1); background: linear-gradient(to bottom, rgba(252,252,252,1) 0%, rgba(237,237,237,1) 100%); border-radius: 72px; box-shadow: 0px 0px 10px rgba(0,0,0,0.2);}"

pyvis :: Graph -> Html ()
pyvis (nodes, edges) = do
  html_ $ do
    head_ $ do
      meta_ [charset_ "utf-8"]
      -- Scripts and links
      script_ [src_ "lib/bindings/utils.js"] ("" :: Text)
      link_
        [ rel_ "stylesheet"
        , href_
            "https://cdnjs.cloudflare.com/ajax/libs/vis-network/9.1.2/dist/dist/vis-network.min.css"
        , integrity_
            "sha512-WgxfT5LWjfszlPHXRmBWHkV2eceiWTOBvrKCNbdgDYTHrT2AeLCGbF4sZlZw3UMN3WtL0tGUoIAKsu8mllg/XA=="
        , crossorigin_ "anonymous"
        , referrerpolicy_ "no-referrer"
        ]
      script_
        [ src_
            "https://cdnjs.cloudflare.com/ajax/libs/vis-network/9.1.2/dist/vis-network.min.js"
        , integrity_
            "sha512-LnvoEWDFrqGHlHmDD2101OrLcbsfkrzoSpvtSQtxK3RMnRV0eOkhhBN2dXHKRrUU8p2DGRTk35n4O8nWSVe1mQ=="
        , crossorigin_ "anonymous"
        , referrerpolicy_ "no-referrer"
        ]
        ("" :: Text)
      link_
        [ rel_ "stylesheet"
        , href_
            "https://cdn.jsdelivr.net/npm/bootstrap@5.0.0-beta3/dist/css/bootstrap.min.css"
        , integrity_
            "sha384-eOJMYsd53ii+scO/bJGFsiCZc+5NDVN2yr8+0RDqr0Ql0h+rP48ckxlpbzKgwra6"
        , crossorigin_ "anonymous"
        ]
      script_
        [ src_
            "https://cdn.jsdelivr.net/npm/bootstrap@5.0.0-beta3/dist/js/bootstrap.bundle.min.js"
        , integrity_
            "sha384-JEW9xMcG8R+pH31jmWH6WWP0WintQrMb4s7ZOdauHnUtxwoG2vI5DkLtS3qm9Ekf"
        , crossorigin_ "anonymous"
        ]
        ("" :: Text)
      -- Inline styles
      style_ [type_ "text/css"] cssStyle
      title_ "Metabolic expansion"
  body_ $ do
    div_ [class_ "card", style_ "width: 100%"]
      $ div_ [id_ "mynetwork", class_ "card-body"] ""
    div_ [id_ "loadingBar"]
      $ div_ [class_ "outerBorder"]
      $ do
          div_ [id_ "text"] "0%"
          div_ [id_ "border"] $ div_ [id_ "bar"] ""
    script_ [type_ "text/javascript"]
      $ toHtmlRaw
          (Data.Text.unlines
             [ "var edges; var nodes; var allNodes; var allEdges; var nodeColors; var originalNodes; var network; var container; var options, data; var filter = { item : '', property : '', value : [] }; function drawGraph() { var container = document.getElementById('mynetwork'); "
             , "nodes = new vis.DataSet("
             , TE.decodeUtf8 . BL.toStrict $ encode nodes
             , ");"
             , "edges = new vis.DataSet("
             , TE.decodeUtf8 . BL.toStrict $ encode edges
             , "); nodeColors = {}; allNodes = nodes.get({ returnType: 'Object' }); for (nodeId in allNodes) { nodeColors[nodeId] = allNodes[nodeId].color; } allEdges = edges.get({ returnType: 'Object' }); data = {nodes: nodes, edges: edges}; var options = { 'configure': { 'enabled': false }, 'edges': { 'color': { 'inherit': true }, 'smooth': { 'enabled': true, 'type': 'dynamic' } }, 'interaction': { 'dragNodes': true, 'hideEdgesOnDrag': false, 'hideNodesOnDrag': false }, 'physics': { 'enabled': true, 'stabilization': { 'enabled': true, 'fit': true, 'iterations': 1000, 'onlyDynamicEdges': false, 'updateInterval': 50 } } }; network = new vis.Network(container, data, options); network.on('stabilizationProgress', function(params) { document.getElementById('loadingBar').removeAttribute('style'); var maxWidth = 496; var minWidth = 20; var widthFactor = params.iterations/params.total; var width = Math.max(minWidth,maxWidth * widthFactor); document.getElementById('bar').style.width = width + 'px'; document.getElementById('text').innerHTML = Math.round(widthFactor*100) + '%'; }); network.once('stabilizationIterationsDone', function() { document.getElementById('text').innerHTML = '100%'; document.getElementById('bar').style.width = '496px'; document.getElementById('loadingBar').style.opacity = 0; setTimeout(function () {document.getElementById('loadingBar').style.display = 'none';}, 500); }); return network; } drawGraph();"
             ])
