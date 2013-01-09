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
package ubic.BAMSandAllen.AllenDataLoaders.Brainspan;

import java.io.FileReader;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import cern.colt.list.DoubleArrayList;

import ubic.BAMSandAllen.Util;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.DescriptiveWithMissing;
import au.com.bytecode.opencsv.CSVReader;

public class BrainSpanPrenatalLMD extends BrainSpanGeneric {
    // protected final String GENE_SYMBOL_KEY = "gene-symbol";

    private static Log log = LogFactory.getLog( BrainSpanPrenatalLMD.class.getName() );

    public BrainSpanPrenatalLMD() throws Exception {
        super();

        GENE_SYMBOL_KEY = "gene-symbol";
        AGE_KEY = "donor_age";

        this.type = "lmacsd1";
        baseFolder = "/home/leon/Downloads/prenatalLMDmicroarray/lmacsd1/";

        // read in columns
        List<Map<String, String>> colMaps = getLines( baseFolder + "Columns.csv" );

        int lineNumber = 0;
        for ( Map<String, String> colMap : colMaps ) {
            colMap.put( "line", lineNumber++ + "" );

        }

        // read in rows
        List<Map<String, String>> rowMaps = getLines( baseFolder + "Rows.csv" );

        log.info( "rowMaps.size(), colMaps.size() " + rowMaps.size() + "x" + colMaps.size() );

        data = new DenseDoubleMatrix<Map<String, String>, Map<String, String>>( rowMaps.size(), colMaps.size() );

        data.setColumnNames( colMaps );

        log.info( "setting row maps" );
        data.setRowNames( rowMaps );
        log.info( "done rownames" );

        String dataFile = baseFolder + "Expression.csv";
        CSVReader reader = new CSVReader( new FileReader( dataFile ) );
        String line[];
        int row = 0;
        while ( ( line = reader.readNext() ) != null ) {
            // log.info( line.length );
            if ( row % 10000 == 0 ) {
                log.info( "Reading row:" + row );
                // break;
            }
            // shift over by one because the first column is the gene ID (probably should check that against the row
            // data
            // for blueprint the first col is row number
            for ( int i = 1; i < line.length; i++ ) {
                data.set( row, i - 1, Double.parseDouble( line[i] ) );
            }
            row++;
        }
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // TODO Auto-generated method stub
        BrainSpanPrenatalLMD lmd = new BrainSpanPrenatalLMD();
        DoubleMatrix<String, String> result = lmd.getRegionByAge( "lmacsd1" );

        String filename = "/home/leon/Downloads/prenatalLMDmicroarray/lmacsd1/matrix.csv";
        Util.writeRTable( filename, result );

        

    }

}
