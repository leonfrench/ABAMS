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
package ubic.BAMSandAllen.gene2region;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.RankedGeneListLoader;
import ubic.BAMSandAllen.Util;
import ubic.basecode.dataStructure.CountingMap;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.io.reader.DoubleMatrixReader;
import ubic.basecode.util.FileTools;

public class GeneToRegionMatrixReader {
    private static Log log = LogFactory.getLog( GeneToRegionMatrixReader.class.getName() );

    String filename;
    DoubleMatrix<String, String> matrix;
    List<GeneToRegion> list;

    public GeneToRegionMatrixReader( String filename ) throws Exception {
        this.filename = filename;
        DoubleMatrixReader reader = new DoubleMatrixReader();
//        reader.setTopLeft( false );
        matrix = reader.read( filename );
        log.info( matrix.rows() + " X " + matrix.columns() );
        List<String> newNames = new LinkedList<String>();
        for ( String rowName : matrix.getRowNames() ) {
            newNames.add( rowName.replaceAll( "\"", "" ) );
        }
        matrix.setRowNames( newNames );
        createList();
    }

    public void createList() {
        list = new LinkedList<GeneToRegion>();
        int nanCount = 0;
        for ( String region : matrix.getColNames() ) {
            for ( String gene : matrix.getRowNames() ) {
                double correlation = matrix.getByKeys( gene, region );
                if ( !Double.isNaN( correlation ) ) {
                    GeneToRegion g2r = new GeneToRegion( gene, region, correlation );
                    list.add( g2r );
                } else {
                    nanCount++;
                }
            }
        }
        Collections.sort( list );
        log.info( "Done sorting list, NaN's = " + nanCount );
    }

    public void printTails( int amountToPrint ) throws Exception {
        log.info( "Printing " + amountToPrint + " from each tail" );
        CountingMap<String> topGenes = new CountingMap<String>();
        CountingMap<String> topRegions = new CountingMap<String>();

        for ( int i = 0; i < amountToPrint; i++ ) {
            GeneToRegion g2r = list.get( i );
            topGenes.increment( g2r.gene );
            topRegions.increment( g2r.region );
            log.info( g2r.toString() );
        }
        for ( int i = 0; i < amountToPrint; i++ ) {
            GeneToRegion g2r = list.get( ( list.size() - 1 ) - i );
            topGenes.increment( g2r.gene );
            topRegions.increment( g2r.region );
            log.info( g2r.toString() );
        }

        String topGenesFile = filename + "." + topGenes.size() + ".topGenes";
        FileTools.stringsToFile( topGenes.sortedKeyList(), topGenesFile );

        RankedGeneListLoader tempLoader = new RankedGeneListLoader( topGenesFile );
        tempLoader.forErmineJ();

        log.info( "________________________" );
        printCountingMap( topGenes );
        log.info( "________________________" );
        printCountingMap( topRegions );

        // log.info( arg0 )
    }

    private void printCountingMap( CountingMap<String> topGenes ) {
        for ( String key : topGenes.sortedKeyList() ) {
            log.info( key + ":" + topGenes.get( key ) );
        }
    }

    public void printTails( double fraction ) throws Exception {
        int amountToPrint = ( int ) ( fraction * list.size() / 2 );
        printTails( amountToPrint );
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        GeneToRegionMatrixReader matrix = new GeneToRegionMatrixReader(
                "/grp/java/workspace/BAMSandAllen/data/ConnectivityAndAllenExpressionMatrixPair.GeneToRegion.jacknife.txt" );
        matrix.printTails( 0.0001 );
    }

}
