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

import java.util.Iterator;

import org.apache.commons.configuration.Configuration;
import org.apache.commons.configuration.PropertiesConfiguration;

public class SetupParameters {

    static String filename = "ABAMS.properties";

    // just grab the config and call get to get the parameter value
    public static Configuration config;

    static {
        try {
            config = new PropertiesConfiguration( filename );
        } catch ( Exception e ) {
            System.out.println( "Could not load " + filename );
            System.exit( 1 );
        }
    }

    public static String getDataFolder() {
        return SetupParameters.config.getString( "abams.dataFolder" );
    }

    public static void main( String argsp[] ) {
        Iterator i = config.getKeys( "abams" );

        while ( i.hasNext() ) {
            String key = ( String ) i.next();
            System.out.println( key + " = " + config.getString( key ) );
        }
    }
}
