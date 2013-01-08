/*******************************************************************************
 * The ABAMS project
 * 
 * Copyright (c) 2012 University of British Columbia
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *       http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package ubic.BAMSandAllen.visualization;

import giny.view.NodeView;

import java.io.FileWriter;
import java.util.Map;

import javax.vecmath.Point3d;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import cytoscape.CyEdge;
import cytoscape.CyNetwork;
import cytoscape.CyNode;
import cytoscape.Cytoscape;
import cytoscape.data.Semantics;
import cytoscape.data.readers.VisualStyleBuilder;
import cytoscape.data.writers.XGMMLWriter;
import cytoscape.view.CyNetworkView;
import cytoscape.visual.VisualPropertyType;

public class FirstTry {
    private static Log log = LogFactory.getLog( FirstTry.class.getName() );
    public static final String FULL_NAME = "fullName";

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // iterate allen regions
        CyNetwork cyNetwork = Cytoscape.createNetwork( "ABAMS", true );

        CyNetworkView view = Cytoscape.createNetworkView( cyNetwork, "network1" );
        // view = Cytoscape.getCurrentNetworkView() ;
        // Cytoscape.getCurrentNetworkView().updateView();
        StructureCatalogLoader oldLoader = new StructureCatalogLoader();
        AllenAtlasAnnotationLoader loader = new AllenAtlasAnnotationLoader();
        Map<String, Point3d> centers = loader.getCenters();
        for ( String region : centers.keySet() ) {
            Point3d point = centers.get( region );
            // log.info( region + ":" + centers.get( region ) );

            String acro = oldLoader.getAcro( region );
            if ( acro == null ) {
                log.info( "Bad acro:" + region );
                System.exit( 1 );
            }
            CyNode node = Cytoscape.getCyNode( acro, true );

            cyNetwork.addNode( node );

            Cytoscape.getNodeAttributes().setAttribute( node.getIdentifier(), FULL_NAME, region );

            NodeView nv = view.addNodeView( node.getRootGraphIndex() );

            nv.setXPosition( point.x * 10 );
            nv.setYPosition( point.y * 10 );

            // NodeView nv = view.getNodeView(node0);
            // view.setNodeIntProperty( node3.getRootGraphIndex(), , 3 );

        }

        // VisualStyleBuilder graphStyle = new VisualStyleBuilder("x", false);
        // graphStyle.setNodeSizeLocked(false);

        // graphStyle.addProperty(node1.getIdentifier(), VisualPropertyType.NODE_LINETYPE, "30");


        CyNode node0 = Cytoscape.getCyNode( "rain", true );
        CyNode node1 = Cytoscape.getCyNode( "rainbow", true );
        CyNode node2 = Cytoscape.getCyNode( "rabbit", true );
        CyNode node3 = Cytoscape.getCyNode( "yellow", true );

        cyNetwork.addNode( node0 );
        cyNetwork.addNode( node1 );
        cyNetwork.addNode( node2 );
        cyNetwork.addNode( node3 );

        CyEdge edge0 = Cytoscape.getCyEdge( node0, node1, Semantics.INTERACTION, "pp", true, true );
        CyEdge edge1 = Cytoscape.getCyEdge( node0, node2, Semantics.INTERACTION, "pp", true, true );
        CyEdge edge3 = Cytoscape.getCyEdge( node2, node1, Semantics.INTERACTION, "pp", true, true );
        cyNetwork.addEdge( edge0 );
        cyNetwork.addEdge( edge1 );
        cyNetwork.addEdge( edge3 );

//        VisualStyleBuilder graphStyle = new VisualStyleBuilder( "x", false );
//        // graphStyle.setNodeSizeLocked(false);
//        graphStyle.addProperty( edge3.getIdentifier(), "weight",  22 );
//        graphStyle.addProperty( edge1.getIdentifier(), VisualPropertyType.EDGE_LINE_WIDTH, "" + 2 );

        Cytoscape.getEdgeAttributes().setAttribute( edge3.getIdentifier(), "w", 2);
        Cytoscape.getEdgeAttributes().setAttribute( edge1.getIdentifier(), "w", 33);
        Cytoscape.getEdgeAttributes().setAttribute( edge0.getIdentifier(), "w", 5);

        log.info( "Edges:" + cyNetwork.getEdgeCount() );
        
        XGMMLWriter writer = new XGMMLWriter( cyNetwork, view );

        FileWriter fw = new FileWriter( SetupParameters.config.getString( "abams.visualization.output" )
                + "FirstTryTemp.xgmml" );
        writer.write( fw );
        fw.close();
        log.info( "Written" );


    }
}
