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

import java.io.FileReader;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.SetupParameters;

import au.com.bytecode.opencsv.CSVReader;

public class BrainStructuresCSVLoader {
    private static Log log = LogFactory.getLog( BrainStructuresCSVLoader.class.getName() );

    Map<String, String> IDtoName;

    public BrainStructuresCSVLoader() throws Exception {
        IDtoName = new HashMap<String, String>();
        CSVReader reader = new CSVReader( new FileReader( SetupParameters.config.getString( "abams.allenAtlasCSV" ) ) );
        String[] line;
        while ( null != ( line = reader.readNext() ) ) {
            if ( line[0].startsWith( "StructureName" ) ) {
                continue;
            }
            String name = line[0].trim();
            String id = line[6].trim();
            //log.info( name + "->" + id );
            IDtoName.put( id, name );
        }
    }

    public String getName( String id ) {
        return IDtoName.get( id );
    }

    public Collection<String> getNames() {
        return IDtoName.values();
    }

    public String getID( String name ) {
        for ( String id : getIDs() ) {
            if ( getName( id ).equals( name ) ) return id;
        }
        return null;
    }

    public Set<String> getIDs() {
        return IDtoName.keySet();
    }

    /**
     * quick and dirty
     * 
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        BrainStructuresCSVLoader loader = new BrainStructuresCSVLoader();
        StructureCatalogLoader oldLoader = new StructureCatalogLoader();

        for ( String name : oldLoader.getRegions() ) {
            String id = loader.getID( name );
            log.info( id );
            if ( id == null )
                log.info( "Null:" + name );
            else
                log.info( loader.getID( name ) );

        }

    }
}
