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

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import ubic.BAMSandAllen.AllenDataLoaders.AllenTop50DataLoader;

import com.sun.org.apache.xalan.internal.xsltc.trax.TrAXFilter;

public class CountTopGenes {

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // TODO Auto-generated method stub
        
        
        AllenTop50DataLoader loader = new AllenTop50DataLoader();

        Map<String, Integer> geneCounts = new HashMap<String, Integer>();

        for ( NomenClatureEntry entry : loader.getAllenEnrichedEntries() ) {
            Set<String> genes = entry.expressedGenes;

            for ( String gene : genes ) {
                int value = 0;
                if ( geneCounts.get( gene ) != null ) value = geneCounts.get( gene );
                geneCounts.put( gene, ++value );
            }
        }

        for ( String gene : geneCounts.keySet() ) {
            System.out.println( gene + "," + geneCounts.get( gene ) );
        }

    }
}
