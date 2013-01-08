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

import java.util.HashSet;
import java.util.Set;

import ubic.BAMSandAllen.AllenDataLoaders.AllenMajorMatrices;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class AllenMatrixExplore {

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {

        Set<String> haveDensityZeros = new HashSet<String>();
        Set<String> haveLevelZeros = new HashSet<String>();
        AllenMajorMatrices allenData = new AllenMajorMatrices();
        

        DoubleMatrix<String, String> levels = allenData.getLevels();
        DoubleMatrix<String, String> density = allenData.getDensities();

        // (DenseDoubleMatrix2DNamed)
        /*
         * FileWriter f = new FileWriter("data\\levelToString.txt"); f.write( levels.toString() ); f.close(); f = new
         * FileWriter("data\\densityToString.txt"); f.write( density.toString() ); f.close();
         */

        double highest = 0;
        for ( int i = 0; i < levels.rows(); i++ ) {
            double sumL = 0;
            double sumD = 0;
            for ( int j = 0; j < levels.columns(); j++ ) {
                sumL += levels.get( i, j );
                sumD += density.get( i, j );
                if ( density.get( i, j ) > highest ) {
                    highest = density.get( i, j );
                }
            }
            if ( sumL == 0 ) {
                try {
                    haveLevelZeros.add( ( String ) levels.getRowName( i ) );
                } catch ( Exception e ) {
                    break;
                }
                // System.out.println(levels.getRowName( i ));
            }
            if ( sumD == 0 ) {
                haveDensityZeros.add( ( String ) levels.getRowName( i ) );
                // System.out.println(levels.getRowName( i ));
            }
        }
        System.out.println( haveLevelZeros.size() + " have all zeroes for levels" );
        System.out.println( haveDensityZeros.size() + " have all zeroes for density" );
        haveLevelZeros.removeAll( haveDensityZeros );
        System.out.println( haveLevelZeros.size() + " are left after removing level zeroes" );

        NeuroNamesMappingLoader NNLoader = new NeuroNamesMappingLoader();
        Set<NomenClatureEntry> swansonEntries = NNLoader.getSwansonEntries();
        Set<NomenClatureEntry> dongEntries = NNLoader.getDongEntries();

        int hit = 0;
        for ( Object s : levels.getColNames() ) {
            String colName = ( String ) s;
            for ( NomenClatureEntry dongEntry : dongEntries ) {
                if ( colName.equalsIgnoreCase( dongEntry.name ) ) {
                    System.out.print( "" + colName + "\t" + dongEntry.NNName + "\tNNID" + dongEntry.NNID );
                    hit++;
                }
            }
            System.out.println( "\t" + colName );
        }
        System.out.println( "hits:" + hit );

        System.out.println( highest );

    }
}
