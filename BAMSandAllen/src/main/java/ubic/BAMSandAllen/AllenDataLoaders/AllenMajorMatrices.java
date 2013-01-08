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
package ubic.BAMSandAllen.AllenDataLoaders;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.xerces.parsers.DOMParser;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

import ubic.BAMSandAllen.SetupParameters;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

/**
 * Generate a level and density matrix from the Allen data
 * 
 * @author lfrench
 */
public class AllenMajorMatrices {

    private DoubleMatrix<String, String> levels;
    private DoubleMatrix<String, String> densities;

    // number of unique ID's minus the number of no gene expression info
    private static final int totalGenes = 20324 - 157;
    private static final int brainRegions = 17;
    public static final String densityFilename = SetupParameters.config.getString( "abams.dataFolder" )
            + "AllenDensities";
    public static final String levelFilename = SetupParameters.config.getString( "abams.dataFolder" ) + "AllenLevels";

    /**
     * @param args
     */
    public AllenMajorMatrices() {
        this( true );
    }

    public AllenMajorMatrices( boolean load ) {
        levels = new DenseDoubleMatrix<String, String>( totalGenes, brainRegions );
        densities = new DenseDoubleMatrix<String, String>( totalGenes, brainRegions );

        try {
            if ( load ) readIn();
        } catch ( Exception e ) {
            System.out.println( "Error loading matrices, please generate them using MakeAllenMatrices.main" );
            e.printStackTrace();
            System.exit( 1 );
        }

    }

    public void writeOut() throws IOException {
        ObjectOutputStream out = new ObjectOutputStream( new FileOutputStream( densityFilename ) );
        out.writeObject( densities );
        out.close();
        out = new ObjectOutputStream( new FileOutputStream( levelFilename ) );
        out.writeObject( levels );
        out.close();
    }

    public void readIn() throws Exception {
        ObjectInputStream in = new ObjectInputStream( new FileInputStream( densityFilename ) );
        densities = ( DenseDoubleMatrix<String, String> ) in.readObject();
        in.close();

        in = new ObjectInputStream( new FileInputStream( levelFilename ) );
        levels = ( DenseDoubleMatrix<String, String> ) in.readObject();
        in.close();
    }

    public void addMeasurement( String geneid, String brainregion, String density, String level ) {
        // check if we have the column/row
        if ( levels.getRowNames().contains( geneid ) ) {
        } else {
            levels.addRowName( geneid );
            densities.addRowName( geneid );
        }

        if ( !levels.getColNames().contains( brainregion ) ) {
            levels.addColumnName( brainregion );
            densities.addColumnName( brainregion );
        }

        levels.setByKeys( geneid, brainregion, Double.parseDouble( level ) );
        densities.setByKeys( geneid, brainregion, Double.parseDouble( density ) );
    }

    public static void main( String[] args ) throws Exception {
        AllenMajorMatrices matrices = new AllenMajorMatrices( false );
        Set<String> seenIDs = new HashSet<String>();

        // get statistics
        File xmlFolder = new File( SetupParameters.config.getString( "abams.dataFolder" ) + "AllenXML/all/" );
        File[] xmlFiles = xmlFolder.listFiles( new FilenameFilter() {
            public boolean accept( File dir, String n ) {
                return n.toLowerCase().endsWith( "xml" );
            }
        } );
        int noGeneExp = 0;
        int noEntrezID = 0;
        int genes = 0;

        for ( File xmlFile : xmlFiles ) {
            InputSource src = new InputSource( new FileInputStream( xmlFile ) );
            DOMParser prsr = new org.apache.xerces.parsers.DOMParser();
            prsr.parse( src );
            Document doc = prsr.getDocument();
            Element element = doc.getDocumentElement();
            System.out.println( xmlFile.getName() + ":" + element.getElementsByTagName( "gene" ).getLength() );
            for ( Element gene : getElementList( element.getElementsByTagName( "gene" ) ) ) {
                genes++;
                Element entrezgeneid = getSingleElement( gene.getElementsByTagName( "entrezgeneid" ) );
                // System.out.println( "Entrez:" + unmarshallText( entrezgeneid ) );
                if ( unmarshallText( entrezgeneid ).equals( "0" ) ) {
                    noEntrezID++;
                }

                Element geneidElement = ( Element ) gene.getElementsByTagName( "geneid" ).item( 0 );
                // System.out.println( "geneid:" + unmarshallText( geneidElement ) );
                if ( seenIDs.contains( unmarshallText( geneidElement ) ) ) {
                    // System.out.println("Seen this one twice:" + unmarshallText( geneidElement ));
                }
                seenIDs.add( unmarshallText( geneidElement ) );

                // should be only one genexpressions tag
                if ( gene.getElementsByTagName( "gene-expressions" ).getLength() == 0 ) {
                    noGeneExp++;
                    continue;
                }
                Element geneExpressions = getSingleElement( gene.getElementsByTagName( "gene-expressions" ) );
                for ( Element geneExpression : getElementList( geneExpressions.getElementsByTagName( "gene-expression" ) ) ) {
                    String avgdensity = null, avglevel = null, structurename = null, geneid = null;
                    for ( Node n : getNodeList( geneExpression.getChildNodes() ) ) {

                        if ( n.getNodeName().equals( "avgdensity" ) ) {
                            avgdensity = unmarshallText( n );
                        } else if ( n.getNodeName().equals( "structurename" ) ) {
                            structurename = unmarshallText( n );
                        } else if ( n.getNodeName().equals( "geneid" ) ) {
                            geneid = unmarshallText( n );
                        } else if ( n.getNodeName().equals( "avglevel" ) ) {
                            avglevel = unmarshallText( n );
                        }
                        // System.out.println( n.getNodeName() + "=" + unmarshallText( n ) );
                    }
                    matrices.addMeasurement( geneid, structurename, avgdensity, avglevel );
                }
            }
        }
        System.out.println( genes + " total gene entries" );
        System.out.println( seenIDs.size() + " gene identifiers" );
        System.out.println( noEntrezID + " have no Entrezgene ID information" );
        System.out.println( noGeneExp + " have no gene expression information" );

        matrices.writeOut();

        // System.out.println(matrices.getDensities().getRowNames());
        // matrices.readIn();
        // System.out.println(matrices.getDensities().toString());

    }

    public static Element getSingleElement( NodeList nodeList ) throws Exception {
        if ( nodeList.getLength() == 1 ) {
            return ( Element ) nodeList.item( 0 );
        } else {
            throw new Exception( "Error node list has more than one or zero" );
        }
    }

    // convience method, make it java 1.5 friendly
    public static List<Element> getElementList( NodeList nodeList ) {
        List<Element> result = new LinkedList<Element>();
        for ( int i = 0; i < nodeList.getLength(); i++ ) {
            result.add( ( Element ) nodeList.item( i ) );
        }
        return result;
    }

    public static List<Node> getNodeList( NodeList nodeList ) {
        List<Node> result = new LinkedList<Node>();
        for ( int i = 0; i < nodeList.getLength(); i++ ) {
            result.add( nodeList.item( i ) );
        }
        return result;
    }

    public static String unmarshallText( Node textNode ) {
        StringBuffer buf = new StringBuffer();

        Node n;
        NodeList nodes = textNode.getChildNodes();

        for ( int i = 0; i < nodes.getLength(); i++ ) {
            n = nodes.item( i );

            if ( n.getNodeType() == Node.TEXT_NODE ) {
                buf.append( n.getNodeValue() );
            } else {
                // expected a text-only node!
            }
        }
        return buf.toString();
    }

    public DoubleMatrix<String, String> getDensities() {
        return densities;
    }

    public DoubleMatrix<String, String> getLevels() {
        return levels;
    }
}
