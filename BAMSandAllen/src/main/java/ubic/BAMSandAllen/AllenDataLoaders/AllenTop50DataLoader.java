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

import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;


import ubic.BAMSandAllen.NeuroNamesMappingLoader;
import ubic.BAMSandAllen.NomenClatureEntry;
import ubic.BAMSandAllen.SetupParameters;
import ubic.basecode.io.excel.ExcelUtil;

public class AllenTop50DataLoader {
    private final String spreadSheetLocation = SetupParameters.config.getString( "abams.dataFolder" )
            + "EnrichedGenesMasterList.xls";

    private Set<String> allenRegions;
    private Set<NomenClatureEntry> allenEnrichedEntries;

    public AllenTop50DataLoader() {
        try {
            NeuroNamesMappingLoader NNLoader = new NeuroNamesMappingLoader();

            allenRegions = new HashSet<String>();

            POIFSFileSystem fs = new POIFSFileSystem( new FileInputStream( spreadSheetLocation ) );
            HSSFWorkbook allen = new HSSFWorkbook( fs );
            HSSFSheet sheet;
            for ( int i = 0; i < allen.getNumberOfSheets(); i++ ) {
                sheet = allen.getSheetAt( i );
                allenRegions.add( allen.getSheetName( i ) );
                // System.out.println( allen.getSheetName( i ) );
            }

            Set<String> allenRegions = getAllenRegions();
            Set<NomenClatureEntry> dongEntries = NNLoader.getDongEntries();
            allenEnrichedEntries = new HashSet<NomenClatureEntry>();

            for ( String allenEnrichedAcro : allenRegions ) {
                for ( NomenClatureEntry dongEntry : dongEntries ) {
                    if ( allenEnrichedAcro.equals( dongEntry.acro ) ) {
                        allenEnrichedEntries.add( dongEntry );
                    }
                }
            }

            // make entry <-> gene links
            for ( NomenClatureEntry dongEntry : allenEnrichedEntries ) {
                sheet = allen.getSheet( dongEntry.acro );
                // System.out.println(dongEntry.acro);
                for ( short i = 2; true; i++ ) {
                    String gene = ExcelUtil.getValue( sheet, i, ( short ) 0 );
                    if ( gene == null ) break;
                    dongEntry.expressedGenes.add( gene );
                    // System.out.println( gene );

                }
            }

        } catch ( Exception e ) {
            e.printStackTrace();
            System.exit( 1 );
        }

    }

    public Set<String> getAllenRegions() {
        return allenRegions;
    }

    public Set<NomenClatureEntry> getAllenEnrichedEntries() throws IOException {
        return allenEnrichedEntries;
    }

    public Set<String> getGenesForDongRegion( String fullName ) {
        for ( NomenClatureEntry dongEntry : allenEnrichedEntries ) {
            if ( dongEntry.name.equalsIgnoreCase( fullName ) ) return dongEntry.expressedGenes;
        }
        return null;
    }

    public static void main( String args[] ) {
        new AllenTop50DataLoader();
    }
}
