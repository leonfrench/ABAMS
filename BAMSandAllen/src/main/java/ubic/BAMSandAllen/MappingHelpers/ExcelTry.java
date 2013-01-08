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
package ubic.BAMSandAllen.MappingHelpers;

import java.io.FileInputStream;
import java.util.HashSet;
import java.util.Set;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;

import ubic.BAMSandAllen.NeuroNamesMappingLoader;
import ubic.BAMSandAllen.NomenClatureEntry;
import ubic.BAMSandAllen.AllenDataLoaders.AllenTop50DataLoader;

public class ExcelTry {

    public static void main( String args[] ) throws Exception {
        // Map<String,String>

        // open neuronames
        NeuroNamesMappingLoader NNLoader = new NeuroNamesMappingLoader();

        Set<NomenClatureEntry> swansonEntries = NNLoader.getSwansonEntries();
        Set<NomenClatureEntry> dongEntries = NNLoader.getDongEntries();
        AllenTop50DataLoader allen = new AllenTop50DataLoader();
        // open Allen enriched genes
        Set<String> allenRegions = allen.getAllenRegions();

        String NNID;
        for ( String allenEnrichedAcro : allenRegions ) {
            // find its entry in Dong
            int hits = 0;
            NomenClatureEntry dongMatch = null;
            for ( NomenClatureEntry dongEntry : dongEntries ) {
                if ( allenEnrichedAcro.equals( dongEntry.acro ) ) {
                    dongMatch = dongEntry;
                    hits++;
                }
            }
            if ( hits != 1 ) {
                System.out
                        .println( "too few or too many matches for a Allen acro:" + hits + " for " + allenEnrichedAcro );
            }
            NNID = dongMatch.NNID;

            // go through Swanson entries and find the NNID linkage
            Set<NomenClatureEntry> matchedSwansonEntries = new HashSet<NomenClatureEntry>();
            for ( NomenClatureEntry swansonEntry : swansonEntries ) {
                if ( NNID.equals( swansonEntry.NNID ) ) {
                    matchedSwansonEntries.add( swansonEntry );
                }
            }
            if ( matchedSwansonEntries.size() == 0 ) {
                System.out.print( "!!" + "NO Swanson MATCHES FOR " );
            }
            System.out.print( "Dong:*" + dongMatch.name + "* maps to Neuronames:" + dongMatch.NNName + "(NNID"
                    + dongMatch.NNID + ")" );
            // gota do multiline, if its just one then stay on the line
            if ( matchedSwansonEntries.size() != 1 ) System.out.println();
            for ( NomenClatureEntry match : matchedSwansonEntries ) {
                System.out.println( " in Swanson:*" + match.name + "*" );
            }
        }

    }
}
