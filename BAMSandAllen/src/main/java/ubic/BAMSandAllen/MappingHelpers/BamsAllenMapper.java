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

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.concurrent.CopyOnWriteArraySet;

import ubic.BAMSandAllen.NeuroNamesMappingLoader;
import ubic.BAMSandAllen.NomenClatureEntry;
import ubic.BAMSandAllen.AllenDataLoaders.AllenTop50DataLoader;

public class BamsAllenMapper {

    private static final String mappingFile = "data\\NNtoSwansonMapping.txt";
    // one Dong region can map to multiple Swanson entries
    private static Map<NomenClatureEntry, Set<String>> mapping;

    // map of swanson to genes
    private static Map<String, Set<String>> genes;

    static {
        try {
            mapping = new HashMap<NomenClatureEntry, Set<String>>();
            genes = new HashMap<String, Set<String>>();

            NeuroNamesMappingLoader NNLoader = new NeuroNamesMappingLoader();

            Set<NomenClatureEntry> swansonEntries = NNLoader.getSwansonEntries();
            Set<NomenClatureEntry> dongEntries = NNLoader.getDongEntries();
            AllenTop50DataLoader allen = new AllenTop50DataLoader();
            // open Allen enriched genes
            Set<String> allenRegions = allen.getAllenRegions();

            BufferedReader f = new BufferedReader( new FileReader( mappingFile ) );
            String line;
            StringTokenizer toke;
            while ( ( line = f.readLine() ) != null ) {
                String split[] = line.split( "->" );
                String neuroname = split[0];
                String swansonName = split[1];

                for ( NomenClatureEntry dongEntry : dongEntries ) {
                    if ( neuroname.equals( dongEntry.NNName ) ) {
                        Set<String> current = mapping.get( dongEntry );
                        if ( current == null ) current = new HashSet<String>();
                        current.add( swansonName.toLowerCase() );
                        mapping.put( dongEntry, current );
                        //System.out.println( dongEntry.name + "->" + swansonName );
                    }
                }

            }
            f.close();

            // now get the gene expression
            AllenTop50DataLoader allenLoader = new AllenTop50DataLoader();
            for ( NomenClatureEntry dongEntry : allenLoader.getAllenEnrichedEntries() ) {
                // get the set of swanson names for this dong entry
                for ( String swansonName : changeDongtoSwanson( dongEntry ) ) {
                    // make an entry with its genes
                    genes.put( swansonName, dongEntry.expressedGenes );
                }
            }

        } catch ( Exception e ) {
            e.printStackTrace();
            System.exit( 1 );
        }
    }

    public static Set<String> getAllSwansonNames() {
        Set<String> result = new HashSet<String>();
        for ( Set<String> values : mapping.values() ) {
            result.addAll( values );
        }
        return result;
    }

    public static Set<String> changeDongtoSwanson( NomenClatureEntry dongEntry ) {
        //had some problems with looking things up using keys as objects
        for (NomenClatureEntry x : mapping.keySet()) {
            if (x.name.equals(dongEntry.name)) return mapping.get(x);
        }
        return null;
//        return mapping.get( dongEntry );
    }

    // get genes - thats it, given a BAMS region
    public static Set<String> getGenesBySwansonName( String swansonName ) {
        return genes.get( swansonName );
    }

    public static void main( String args[] ) {
        for ( String test : getAllSwansonNames() ) {
            System.out.print( test );
            System.out.println( getGenesBySwansonName( test ) );
        }
    }
}
