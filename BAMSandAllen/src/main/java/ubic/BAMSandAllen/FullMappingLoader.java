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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;

import ubic.basecode.io.excel.ExcelUtil;



public class FullMappingLoader {
    // colums in the Neuronames Excel file that list acronym and neuronames ID - this is checked
    final static short ALLENCOL = 0; // full name
    final static short BAMSCOL = 4;
    HSSFSheet sheet;
    HSSFWorkbook excelFile;
    Map<String, Set<String>> mapping;

    public FullMappingLoader() {
        try {
            mapping = new HashMap<String, Set<String>>();
            POIFSFileSystem fs2 = new POIFSFileSystem( new FileInputStream( SetupParameters.config
                    .getString( "abams.allNomenclatures" ) ) );
            excelFile = new HSSFWorkbook( fs2 );
            sheet = excelFile.getSheetAt( 0 );
            for ( int i = 1; i < sheet.getLastRowNum() + 1; i++ ) {
                // last row num doesnt seem to work
                if ( ExcelUtil.getValue( sheet, i, ALLENCOL ) == null ) break;
                String dongName = ExcelUtil.getValue( sheet, i, ALLENCOL );
                String BAMSName = ExcelUtil.getValue( sheet, i, BAMSCOL );
                Set<String> current = mapping.get( BAMSName );
                if ( current == null ) current = new HashSet<String>();
                current.add( dongName );
                mapping.put( BAMSName, current );
                // System.out.println(dongName +" -> " + BAMSName);
            }
        } catch ( Exception e ) {
            e.printStackTrace();
            System.exit( 1 );
        }
    }

    public Map<String, Set<String>> getAllentoBamsMapping() {
        return mapping;
    }

}
