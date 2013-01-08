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

import java.io.FileInputStream;
import java.util.HashSet;
import java.util.Set;

import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;

import ubic.basecode.io.excel.ExcelUtil;



public class NeuroNamesMappingLoader {
    // colums in the Neuronames Excel file that list acronym and neuronames ID - this is checked
    final static short NNFULLCOL = 0; // full name
    final static short NNACROCOL = 1;
    final static short NNIDCOL = 4;
    final static short NNDEFAULTNAME = 2;
    HSSFSheet dong;
    HSSFSheet swanson;
    HSSFWorkbook neuronames;
    Set<NomenClatureEntry> dongEntries,swansonEntries;

    public NeuroNamesMappingLoader() {
        try {
            POIFSFileSystem fs2 = new POIFSFileSystem( new FileInputStream( SetupParameters.config.getString( "abams.neuronames.location" ) ) );
            neuronames = new HSSFWorkbook( fs2 );
            dong = neuronames.getSheet( "Dong" );
            swanson = neuronames.getSheet( "Swanson" );
            dongEntries = getEntries( dong );
            swansonEntries = getEntries( swanson );
            checkCols( dong, swanson );
        } catch ( Exception e ) {
            e.printStackTrace();
            System.exit( 1 );
        }
    }


    public static void checkCols( HSSFSheet dong, HSSFSheet swanson ) throws Exception {
        if ( !ExcelUtil.getValue( swanson, 0, NNACROCOL ).equals( "Swn Acro" ) )
            throw new Exception( "Excel columns not matching" );
        ;
        if ( !ExcelUtil.getValue( swanson, 0, NNIDCOL ).equals( "NN ID " ) ) throw new Exception( "Excel columns not matching" );
        ;

        if ( !ExcelUtil.getValue( dong, 0, NNACROCOL ).equals( "Dong Acro" ) )
            throw new Exception( "Excel columns not matching" );
        ;
        if ( !ExcelUtil.getValue( dong, 0, NNIDCOL ).equals( "NN ID " ) ) throw new Exception( "Excel columns not matching" );
    }

    public static Set<NomenClatureEntry> getEntries( HSSFSheet NNsheet ) {
        Set<NomenClatureEntry> result = new HashSet<NomenClatureEntry>();
        NomenClatureEntry e;
        for ( int i = 1; i < NNsheet.getLastRowNum() + 1; i++ ) {
            e = new NomenClatureEntry();
            e.acro = ExcelUtil.getValue( NNsheet, i, NNACROCOL );
            e.name = ExcelUtil.getValue( NNsheet, i, NNFULLCOL );
            e.NNID = ExcelUtil.getValue( NNsheet, i, NNIDCOL );
            e.NNName = ExcelUtil.getValue( NNsheet, i, NNDEFAULTNAME );
            result.add( e );
        }
        return result;
    }

    public Set<NomenClatureEntry> getDongEntries() {
        return dongEntries;
    }

    public Set<NomenClatureEntry> getSwansonEntries() {
        return swansonEntries;
    }
}
