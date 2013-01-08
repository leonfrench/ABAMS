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
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Point3d;

import org.apache.commons.lang.ObjectUtils.Null;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.RankedGeneListLoader;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.MatrixPairs.AllenMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenDataPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.MatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.adjacency.IdentityAdjacency;
import ubic.BAMSandAllen.adjacency.SingleRowAdjacency;
import ubic.BAMSandAllen.adjacency.VoxelVolumeAdjacency;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import cytoscape.CyEdge;
import cytoscape.CyNetwork;
import cytoscape.CyNode;
import cytoscape.Cytoscape;
import cytoscape.data.Semantics;
import cytoscape.data.readers.VisualStyleBuilder;
import cytoscape.data.writers.XGMMLWriter;
import cytoscape.view.CyNetworkView;
import cytoscape.visual.VisualPropertyType;

public class GraphFromABAMS {
    private static Log log = LogFactory.getLog( FirstTry.class.getName() );
    public static final String FULL_NAME = "fullName";
    public static final String DIVISION = "division";
    public static final String EDGE_WEIGHT_NAME = "weight";

    CyNetwork cyNetwork;
    CyNetworkView view;
    String filename;
    StructureCatalogLoader oldLoader;
    AllenAtlasAnnotationLoader loader;

    public GraphFromABAMS() throws Exception {
        cyNetwork = Cytoscape.createNetwork( "ABAMS", true );
        view = Cytoscape.createNetworkView( cyNetwork, "network1" );
        filename = SetupParameters.config.getString( "abams.visualization.output" ) + "Cyto";
        oldLoader = new StructureCatalogLoader();
        loader = new AllenAtlasAnnotationLoader();
    }

    public void addNodes( Collection<String> regionNames ) {
        Map<String, Point3d> centers = loader.getCenters();
        Set<String> divisions = new HashSet<String>();
        divisions.add( "Hindbrain" );
        divisions.add( "Interbrain" );
        divisions.add( "Midbrain" );
        divisions.add( "Cerebrum" );
        divisions.add( "Cerebellum" );

        for ( String region : regionNames ) {
            Point3d point = centers.get( region );
            if ( point == null ) {
                log.info( region + " has null point" );
                continue;
            }
            // log.info( region + ":" + centers.get( region ) );

            String acro = oldLoader.getAcro( region );
            if ( acro == null ) {
                log.info( "Bad acro:" + region );
                System.exit( 1 );
            }

            CyNode node = Cytoscape.getCyNode( acro, true );

            cyNetwork.addNode( node );

            Cytoscape.getNodeAttributes().setAttribute( node.getIdentifier(), FULL_NAME, region );

            for ( String division : divisions ) {
                Set<String> parents = oldLoader.getParents( region );
                // if it is in that division, or is that division then mark it
                if ( parents != null && ( parents.contains( division ) || division.equals( region ) ) ) {
                    Cytoscape.getNodeAttributes().setAttribute( node.getIdentifier(), DIVISION, division );
                }
            }

            NodeView nv = view.addNodeView( node.getRootGraphIndex() );

            nv.setXPosition( point.x * 10 );
            nv.setYPosition( point.y * 10 );

            // NodeView nv = view.getNodeView(node0);
            // view.setNodeIntProperty( node3.getRootGraphIndex(), , 3 );

        }
        log.info( "Node count:" + cyNetwork.getNodeCount() );

    }

    public void writeOut() throws Exception {
        writeOut( "" );
    }

    public void writeOut( String endFix ) throws Exception {
        XGMMLWriter writer = new XGMMLWriter( cyNetwork, view );
        String finalFilename = filename + "." + endFix + ".xgmml";
        FileWriter fw = new FileWriter( finalFilename );
        writer.write( fw );
        fw.close();
        log.info( "Written to " + finalFilename );
    }

    public CyNode getNodeFromAllenRegionName( String name ) {
        return Cytoscape.getCyNode( oldLoader.getAcro( name ) );
    }

    public void deleteNodeFromAllenRegionName( String name ) {
        CyNode node = getNodeFromAllenRegionName( name );
        int node_index = cyNetwork.getIndex( node );
        cyNetwork.removeNode( node_index, true );
    }

    public DoubleMatrix<String, String> thresholdMatrix( DoubleMatrix<String, String> matrix, double fraction ) {
        DoubleMatrix<String, String> result = matrix.copy();
        LinkedList<Double> allValues = new LinkedList<Double>();
        for ( int i = 0; i < matrix.rows(); i++ ) {
            for ( int j = 0; j < matrix.columns(); j++ ) {
                allValues.add( matrix.get( i, j ) );
            }
        }
        Collections.sort( allValues );
        double threshold = allValues.get( ( int ) ( fraction * allValues.size() ) );
        log.info( "threshold:" + threshold );
        for ( int i = 0; i < matrix.rows(); i++ ) {
            for ( int j = 0; j < matrix.columns(); j++ ) {
                if ( result.get( i, j ) < threshold ) {
                    result.set( i, j, 0d );
                }
            }
        }
        return result;
    }

    private void addBAMSMatrix( DoubleMatrix<String, String> matrix, ConnectivityAndAllenDataPair pair, String name,
            boolean symetric ) throws Exception {
        int beforeEdges = cyNetwork.getEdgeCount();
        filename += "." + name;
        cyNetwork.setTitle( cyNetwork.getTitle() + " " + name );
        // add the nodes to the graph
        Set<String> regionNamesBAMS = new HashSet<String>( matrix.getColNames() );
        // add in the rownames too???
        // if it is not using adjacency (like a connection matrix)
        regionNamesBAMS.addAll( matrix.getRowNames() );

        Set<String> regionNamesAllen = new HashSet<String>();
        for ( String BAMSName : regionNamesBAMS ) {
            Set<String> allenMapped = pair.convertANametoB( BAMSName );
            if ( allenMapped.size() > 1 ) {
                log.info( BAMSName + " Has more than one Allen mapping" );
            }
            if ( allenMapped.size() == 0 ) {
                log.info( BAMSName + " Has zero Allen mapping!" );
            }
            regionNamesAllen.addAll( allenMapped );
        }
        addNodes( regionNamesAllen );
        log.info( "Total Allen regions:" + regionNamesAllen.size() );

        DoubleMatrix<String, String> finalMatrix = matrix;

        int calledAdded = 0;
        for ( String iBAMS : finalMatrix.getRowNames() ) { // column
            for ( String iAllen : pair.convertANametoB( iBAMS ) ) {
                for ( String jBAMS : finalMatrix.getColNames() ) { // row
                    for ( String jAllen : pair.convertANametoB( jBAMS ) ) {
                        double value = finalMatrix.getByKeys( iBAMS, jBAMS );
                        int iIndex = finalMatrix.getRowIndexByName( iBAMS );
                        int jIndex = finalMatrix.getColIndexByName( jBAMS );
                        addEdgeFromMatrix( name, symetric, finalMatrix, iAllen, iIndex, jAllen, jIndex, value );
                    }
                }
            }
        }
        log.info( "Added edges:" + ( cyNetwork.getEdgeCount() - beforeEdges ) );

    }

    // example of non-directed is nomenclature
    private void addAllenMatrix( DoubleMatrix<String, String> matrix, String name, boolean symetric ) {
        // add the nodes to the graph
        List<String> regionNames = matrix.getColNames();
        addNodes( regionNames );

        DoubleMatrix<String, String> finalMatrix = matrix;

        // use the adjacency? -yes
        for ( String i : finalMatrix.getRowNames() ) {
            for ( String j : finalMatrix.getColNames() ) {
                double value = finalMatrix.getByKeys( i, j );
                int iIndex = finalMatrix.getRowIndexByName( i );
                int jIndex = finalMatrix.getColIndexByName( j );
                addEdgeFromMatrix( name, symetric, finalMatrix, i, iIndex, j, jIndex, value );
            }
        }
        log.info( "Edges:" + cyNetwork.getEdgeCount() );
    }

    private void addEdgeFromMatrix( String interactionType, boolean symetric, DoubleMatrix<String, String> adjacency,
            String i, int iIndex, String j, int jIndex, double value ) {
        // only make it if it's a nonzero edge value
        if ( value != 0d ) {
            boolean create = true;
            // if it's directed/adjaceny matrix use the whole matrix, undirected -> symetric matrix
            if ( !( symetric && ( iIndex < jIndex ) ) ) {

                CyNode iNode = getNodeFromAllenRegionName( i );
                CyNode jNode = getNodeFromAllenRegionName( j );

                // they probably don't have coordinates and were not made
                if ( iNode == null || jNode == null ) {
                    log.info( i + " -> " + j + " edge not made, null nodes" );
                    return;
                }

                CyEdge edge = Cytoscape.getCyEdge( iNode, jNode, Semantics.INTERACTION, interactionType, create,
                        !symetric );
                cyNetwork.addEdge( edge );
                Cytoscape.getEdgeAttributes().setAttribute( edge.getIdentifier(), EDGE_WEIGHT_NAME, value );

            }
        }
    }

    public void addExpression( ABAMSDataMatrix matrix, String rowName ) {
        int nullNodes = 0;
        for ( String regionName : matrix.getColNames() ) {
            CyNode regionNode = getNodeFromAllenRegionName( regionName );

            // they probably don't have coordinates and were not made
            if ( regionNode == null ) {
                nullNodes++;
                continue;
            }

            Cytoscape.getNodeAttributes().setAttribute( regionNode.getIdentifier(), rowName,
                    matrix.getByKeys( rowName, regionName ) );
        }
        log.info( nullNodes + " null nodes, not made." );
    }

    public void addMatrix( ConnectivityAndAllenExpressionMatrixPair pair, boolean useAdjacency, ABAMSDataMatrix matrix,
            boolean BAMSSpace ) throws Exception {
        boolean threshold = false;
        boolean useSquare = false;
        addMatrix( pair, useAdjacency, useSquare, threshold, 0d, matrix, BAMSSpace );
    }

    public DoubleMatrix<String, String> addMatrix( ConnectivityAndAllenExpressionMatrixPair pair, boolean useAdjacency,
            boolean useSquare, boolean threshold, double thresholdValue, ABAMSDataMatrix matrix, boolean BAMSSpace )
            throws Exception {
        String name;
        DoubleMatrix<String, String> matrixPlain;
        boolean symetric;
        name = matrix.getName();
        // some strange bug in cytoscape when you have brakets
        name = name.replace( '(', '.' );
        name = name.replace( ')', '.' );
        if ( useAdjacency ) {
            name += ".adjacency";
            matrixPlain = matrix.getAdjacency();
            symetric = true;
        } else if ( useSquare ) {
            name += ".square";
            matrixPlain = matrix.getSquare();
            symetric = pair.getDirection().equals( Direction.ANYDIRECTION );
            // symetric = false;
        } else {
            matrixPlain = matrix;
            symetric = false;
        }

        if ( threshold ) {
            name += "." + thresholdValue + ".threshold";
            matrixPlain = thresholdMatrix( matrixPlain, thresholdValue );
        }

        if ( BAMSSpace )
            addBAMSMatrix( matrixPlain, pair, name, symetric );
        else
            addAllenMatrix( matrixPlain, name, symetric );

        return matrixPlain;
    }

    public void addNomenclature() {
        DoubleMatrix<String, String> nomen = oldLoader.getParentMatrix();
        boolean symetric = false;
        addAllenMatrix( nomen, ".Nomenclature", symetric );
    }

    public void addVolume() throws Exception {
        // volume
        AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
        DoubleMatrix<String, String> volMatrix = spaceLoader.getVolumeMatrix();
        ABAMSDataMatrix volABAMS = new ABAMSDataMatrix( volMatrix, "Volume", new VoxelVolumeAdjacency() );
        addExpression( volABAMS, "volume" );
    }

    public void intersectEdges( ConnectivityAndAllenExpressionMatrixPair pair, boolean symetric,
            DoubleMatrix<String, String> matrixAPlain, DoubleMatrix<String, String> matrixBPlain ) throws Exception {
        // make new matrix with names joined

        DoubleMatrix<String, String> matrixHadaProd;
        Set<String> colNames = new HashSet<String>();
        colNames.addAll( matrixAPlain.getColNames() );
        colNames.retainAll( matrixBPlain.getColNames() );
        matrixHadaProd = new DenseDoubleMatrix<String, String>( colNames.size(), colNames.size() );
        matrixHadaProd.setRowNames( new LinkedList<String>( colNames ) );
        matrixHadaProd.setColumnNames( new LinkedList<String>( colNames ) );

        log.info( colNames.size() + " colnames kept of " + matrixAPlain.columns() + " and " + matrixBPlain.columns() );
        int edgeCount = 0;
        for ( String rowName : colNames ) {
            for ( String colName : colNames ) {
                double value = matrixBPlain.getByKeys( rowName, colName ) * matrixAPlain.getByKeys( rowName, colName );
                matrixHadaProd.setByKeys( rowName, colName, value );
                if ( value != 0 ) edgeCount++;
            }
        }
        log.info( "Intersection edges added:" + edgeCount + " (includes doubles)" );

        // DoubleMatrix<String, String> matrixHadaProd = matrixBPlain.copy();
        // if ( !matrixBPlain.getColNames().equals( matrixAPlain.getColNames() ) )
        // throw new Exception( "Error col names don't match" );
        // int edgeCount = 0;
        // for ( int i = 0; i < matrixBPlain.rows(); i++ ) {
        // for ( int j = 0; j < matrixBPlain.columns(); j++ ) {
        // double value = matrixBPlain.get( i, j ) * matrixAPlain.get( i, j );
        // matrixHadaProd.set( i, j, value );
        // if ( value != 0 ) edgeCount++;
        // }
        // }
        // log.info( "Intersection edges added:" + edgeCount + " (includes doubles)" );
        addBAMSMatrix( matrixHadaProd, pair, "intersection", symetric );
    }

    public void setToRegions( Set<String> regionsToSet ) {
        // deleteNodeFromAllenRegionName( name );

    }

    public void addInterstingGenes( Direction direction ) throws Exception {
        boolean removeNonExp = false;
        ABAMSDataMatrix expressionAllenSpace = ExpressionMatrixPairFactory.getEnergyMatrix( direction, removeNonExp );

        GraphFromABAMS graph = this;
        graph.addExpression( expressionAllenSpace, "Pgrmc1[797]" );
        graph.addExpression( expressionAllenSpace, "Rpp25[583571]" );
        graph.addExpression( expressionAllenSpace, "Sema3a[79761021]" );
        graph.addExpression( expressionAllenSpace, "Vamp2[1093]" );
        graph.addExpression( expressionAllenSpace, "Pcp2[77413702]" );
        graph.addExpression( expressionAllenSpace, "Gpc3[71020431]" );
        graph.addExpression( expressionAllenSpace, "Gm47[70565879]" );
        graph.addExpression( expressionAllenSpace, "Nrp2[80514091]" ); // interbrain
        graph.addExpression( expressionAllenSpace, "Kirrel1[71613657]" ); // cerebrum
        graph.addExpression( expressionAllenSpace, "A930033C23Rik*[74300717]" ); // cerebrum -.33, interbrain
        graph.addExpression( expressionAllenSpace, "Tac2[77279001]" ); // interbrain
        graph.addExpression( expressionAllenSpace, "Cpne5[544709]" ); // midbrain 0.14
        graph.addExpression( expressionAllenSpace, "Svil[77332748]" ); // global
        graph.addExpression( expressionAllenSpace, "Ptprr[74882784]" );// global
        graph.addExpression( expressionAllenSpace, "L1cam[80342072]" );// global

        graph.addExpression( expressionAllenSpace, "Nef3[69782174]" );// neurofilament
        graph.addExpression( expressionAllenSpace, "Nef3[73817434]" );// neurofilament
        graph.addExpression( expressionAllenSpace, "Nefh[71920442]" );// neurofilament
        graph.addExpression( expressionAllenSpace, "Nefh[74512048]" );// neurofilament
        graph.addExpression( expressionAllenSpace, "Nefl[70919083]" );// neurofilament
        graph.addExpression( expressionAllenSpace, "Nefl[73512198]" );// neurofilament

        // from addition of genes and direct connectivity fishing
        graph.addExpression( expressionAllenSpace, "Lrrc8d[385997]" );
        graph.addExpression( expressionAllenSpace, "Thsd4[69527776]" );
        graph.addExpression( expressionAllenSpace, "Fhad1[68546125]" );
        graph.addExpression( expressionAllenSpace, "Epha8[69672090]" );
        graph.addExpression( expressionAllenSpace, "Rgs3[1822]" );
        graph.addExpression( expressionAllenSpace, "Gpc3[71020431]" );

        graph.addExpression( expressionAllenSpace, "Hpcal4[73520985]" );
        graph.addExpression( expressionAllenSpace, "Mef2c[79567505]" );

    }

    public static void main( String[] args ) throws Exception {

        Direction direction = AnalyzeBAMSandAllenGenes.Direction.ANYDIRECTION;

        // cannot use virtual regions because it will lead to a non-square connectivity matrix (incoming only) - virtual
        // regions only formed on the columns, not rows
        boolean useVirtual = true;
        boolean removeNonExp = true;

        boolean squareMatrix = true;
        boolean logDistance = true;
        boolean slow = false;
        RegressMatrix regressType = RegressMatrix.BOTH;

        ConnectivityAndAllenExpressionMatrixPair pair = null;

        boolean BAMSSpace = true;
        boolean useAdjacency = false;
        boolean threshold = false;
        boolean useSquare = true;
        boolean run = true;
        double thresholdValue = 0.9;
        ABAMSDataMatrix matrix = null;
        GraphFromABAMS graph = null;
        DoubleMatrix<String, String> matrixAPlain = null;
        graph = new GraphFromABAMS();

        pair = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp, run,
                squareMatrix );
        log.info( "Connections:" + pair.getConnectionCount() );
        pair.printConnectionInfo();
        // set matrix here
        matrix = pair.getMatrixA(); // connectivity
        matrixAPlain = graph.addMatrix( pair, useAdjacency, useSquare, threshold, thresholdValue, matrix, BAMSSpace );

        graph.addInterstingGenes( direction );

        // just write connectivity
        // graph.writeOut( "AnyDirection" );

        // pair.writeRMatrices();
        // pair.writeImages();

        // System.exit( 1 );

        // random connections
        // for ( int i = 1; i < 6; i++ ) {
        // // pair = ExpressionMatrixPairFactory.connectivityPartial( direction, slow, regressType, useVirtual,
        // // removeNonExp,
        // // logDistance );
        //
        // //pair.shuffleConnectivityCols( i );
        //
        // // set matrix here
        // matrix = pair.getMatrixA(); // connectivity
        // matrixAPlain = graph.addMatrix( pair, useAdjacency, threshold, thresholdValue, matrix, BAMSSpace );
        //
        // // just write connectivity
        // graph.writeOut( "Incoming.random.seed" + i );
        // }
        // System.exit( 1 );

        // matrix = pair.getMatrixB();
        // DoubleMatrix<String, String> matrixBPlain = graph.addMatrix( pair, useAdjacency, threshold, thresholdValue,
        // matrix, BAMSSpace );

        direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;

        squareMatrix = false;

        BAMSSpace = true;
        useAdjacency = true;
        threshold = true;
        useSquare = false;
        thresholdValue = 0.9;

        // pair = ExpressionMatrixPairFactory.connectivityPartial( direction, slow, regressType, useVirtual,
        // removeNonExp,
        // logDistance );
        //
        // // get topgenes
        // String filename = "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final ammon/";
        // // filename += "LOOGenesInOrder.in.partialcon.ammon.txt.452.0.01985.topGenes.txt";
        // filename += "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt";
        // // filename += "LOOGenesInOrder.space.ammon.txt.420.0.018435.topGenes.txt";
        // RankedGeneListLoader aLook = new RankedGeneListLoader( filename );
        // pair.setMatrixBDataRows( aLook.getLines() );
        // log.info( pair.getCorrelation() );
        // pair.getMatrixB().setName( "TopOutgoingExpression" );

        // matrix = pair.getMatrixB();
        // useAdjacency = true;
        // useSquare = false;
        // DoubleMatrix<String, String> matrixBPlain = graph.addMatrix( pair, useAdjacency, useSquare, threshold,
        // thresholdValue, matrix, BAMSSpace );
        //
        // matrix = pair.getMatrixA(); // connectivity
        // matrixAPlain = graph.addMatrix( pair, useAdjacency, useSquare, threshold, thresholdValue, matrix, BAMSSpace
        // );
        //
        // graph.intersectEdges( pair, useAdjacency, matrixAPlain, matrixBPlain );
        boolean add = true;
        threshold = false;
        useSquare = true;

        direction = AnalyzeBAMSandAllenGenes.Direction.ANYDIRECTION;
        pair = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp );

        // Mef2c[79567505]
        // Nefh[74512048]
        String rowName = "Hpcal4[73520985]";
        pair.setToSingleGene( rowName, add );

        String operation = "sum";
        if ( !add ) operation = "diff";
        pair.getMatrixB().setName( rowName + "." + operation );

        DoubleMatrix<String, String> singleGene = graph.addMatrix( pair, useAdjacency, useSquare, threshold,
                thresholdValue, pair.getMatrixB(), BAMSSpace );

        // graph.addNomenclature();
        // graph.addVolume();
        boolean undirected = true;
        graph.intersectEdges( pair, undirected, singleGene, matrixAPlain );

        pair.getFlattenedCorrelation( false );

        graph.writeOut( "and.Connections" );
    }
}
