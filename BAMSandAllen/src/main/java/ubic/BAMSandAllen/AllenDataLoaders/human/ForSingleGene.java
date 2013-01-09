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

package ubic.BAMSandAllen.AllenDataLoaders.human;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import cern.colt.list.DoubleArrayList;
import cern.colt.matrix.DoubleMatrix1D;

import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.Brainspan.BrainSpanPrenatalLMD;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.io.reader.DoubleMatrixReader;

public class ForSingleGene {
    private static Log log = LogFactory.getLog( ForSingleGene.class.getName() );

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        String gene = "\"YourGene\"";
        DoubleMatrixReader reader = new DoubleMatrixReader();

        DoubleMatrix<String, String> h9861 = reader
                .read( "/home/leon/Desktop/temp/three human/9861.matrix.29192 x 323.txt" );
        log.info( "Read one" );
        h9861.getRowIndexByName( gene );
        DoubleMatrix<String, String> h10021 = reader
                .read( "/home/leon/Desktop/temp/three human/10021.matrix.29192 x 346.txt" );
        log.info( "Read one" );
        DoubleMatrix<String, String> h12876 = reader
                .read( "/home/leon/Desktop/temp/three human/12876.matrix.29192 x 158.txt" );
        log.info( "Read one" );

        Set<String> regions = new HashSet<String>();
        regions.addAll( h9861.getColNames() );
        regions.addAll( h10021.getColNames() );
        regions.addAll( h12876.getColNames() );

        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( regions.size(), 3 );
        result.setRowNames( new LinkedList<String>( regions ) );
        result.addColumnName( "h9861" );
        result.addColumnName( "h10021" );
        result.addColumnName( "h12876" );

        for ( String region : regions ) {
            log.info( "Region:" + region );

            if ( h9861.getColNames().contains( region ) )
                result.setByKeys( region, "h9861", h9861.getByKeys( gene, region ) );
            else
                result.setByKeys( region, "h9861", Double.NaN );
            if ( h10021.getColNames().contains( region ) )
                result.setByKeys( region, "h10021", h10021.getByKeys( gene, region ) );
            else
                result.setByKeys( region, "h10021", Double.NaN );
            if ( h12876.getColNames().contains( region ) )
                result.setByKeys( region, "h12876", h12876.getByKeys( gene, region ) );
            else
                result.setByKeys( region, "h12876", Double.NaN );
        }

        String filename = "/home/leon/Desktop/temp/three human/human.array.csv";
        Util.writeRTable( filename, result );

    }
}
