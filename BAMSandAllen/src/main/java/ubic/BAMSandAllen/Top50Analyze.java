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
package ubic.BAMSandAllen;

import static ubic.BAMSandAllen.Vocabulary.hasAllenEnrichedGene;

import java.util.List;
import java.util.Set;
import java.util.concurrent.CopyOnWriteArraySet;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ClassSelectors.Top50ClassSelector;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.graphics.MatrixDisplay;

import com.hp.hpl.jena.ontology.OntClass;
import com.hp.hpl.jena.rdf.model.Literal;
import com.hp.hpl.jena.rdf.model.NodeIterator;
import com.hp.hpl.jena.rdf.model.RDFNode;

public class Top50Analyze extends AnalyzeBAMSandAllenGenes {
    private static Log log = LogFactory.getLog( Top50Analyze.class.getName() );

    public Top50Analyze() {
        super( new Top50ClassSelector() );
        // super( new Top50ClassSelectorMinusBST() );
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        Direction direction = Direction.INCOMING;

        Top50Analyze top50 = new Top50Analyze();
        top50.propagateConnections();

        DoubleMatrix<String, String> connectionMatrix;
        connectionMatrix = top50.makeConnectionMatrix( direction );
        // connectionMatrix.makePNG( "data\\connectionMatrixTop50.png", 12 );

        MatrixDisplay matDisplay = new MatrixDisplay( connectionMatrix );
        matDisplay.saveImage( "data\\connectionMatrixTop50.png" );

        DoubleMatrix<String, String> intersectionMatrix = top50.makeIntersectionMatrix();

        // intersectionMatrix.makePNG( "data\\intersectionMatrixTop50.png", 12 );
        matDisplay = new MatrixDisplay( intersectionMatrix );
        matDisplay.saveImage( "data\\intersectionMatrixTop50.png" );

        top50.showMappings();

        for ( int r = 0; r < connectionMatrix.rows(); r++ ) {
            for ( int c = 0; c < connectionMatrix.columns(); c++ ) {
                if ( connectionMatrix.get( r, c ) == 1 ) {
                    System.out.println( "connected," );
                } else {
                    System.out.println( "not connected," );

                }
            }
        }

        // Map<String, Integer> geneCounts;
        // geneCounts = new HashMap<String, Integer>();
        // Double strength = connectionMatrix.getByNames( classToString( regionA ), classToString( regionB ) );
        // if ( strength > 0 ) { //
        // int shared = countSharedGenes( regionA, regionB );
        // System.out.println( "connected|" + shared );
        // } else {
        // System.out.println( "not|" + shared );
        // }
        // System.out.println( "-------------------------------" );
        // for ( String gene : geneCounts.keySet() ) {
        // System.out.println( gene + "," + geneCounts.get( gene ) );
        // }
        // for ( String gene : result ) {
        // int value = 0;
        // if ( geneCounts.get( gene ) != null ) value = geneCounts.get( gene );
        // geneCounts.put( gene, ++value );
        // }

    }

    public Set<String> getTop50Genes( OntClass o ) {
        Set<String> result = new CopyOnWriteArraySet<String>();
        NodeIterator nodes = o.listPropertyValues( hasAllenEnrichedGene );
        for ( Object nodeObj : nodes.toList() ) {
            RDFNode node = ( RDFNode ) nodeObj;
            Literal l = ( Literal ) node.as( Literal.class );
            result.add( l.getString() );
        }
        return result;
    }

    public Set<String> getSharedGenes( OntClass a, OntClass b ) {
        Set<String> result = getTop50Genes( a );
        result.retainAll( getTop50Genes( b ) );
        return result;
    }

    public int countSharedGenes( OntClass a, OntClass b ) {
        return getSharedGenes( a, b ).size();
    }

    public DoubleMatrix<String, String> makeIntersectionMatrix() {
        // only works with the top50
        Set<OntClass> regions = JenaUtil.fitlerClasses( BAMS, new Top50ClassSelector() );
        // get strings
        List<String> enrichedRegionNames = BAMSData.convertRegionstoNames( regions );

        DoubleMatrix<String, String> intersectionMatrix = new DenseDoubleMatrix<String, String>( enrichedRegionNames
                .size(), enrichedRegionNames.size() );
        intersectionMatrix.setRowNames( enrichedRegionNames );
        intersectionMatrix.setColumnNames( enrichedRegionNames );

        for ( OntClass regionA : regions ) {
            for ( OntClass regionB : regions ) {
                if ( regionA.equals( regionB ) ) continue;
                String regionAName = BAMSData.convertClassRegionToString( regionA );
                String regionBName = BAMSData.convertClassRegionToString( regionB );
                // System.out.print( regionAName + "|" + regionBName + "|" );
                int shared = countSharedGenes( regionA, regionB );
                intersectionMatrix.setByKeys( regionBName, regionAName, ( double ) shared );
                intersectionMatrix.setByKeys( regionAName, regionBName, ( double ) shared );
            }
        }
        return intersectionMatrix;
    }

}
