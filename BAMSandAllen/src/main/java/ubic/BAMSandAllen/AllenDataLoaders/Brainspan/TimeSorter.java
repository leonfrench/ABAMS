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

import java.util.Comparator;

import com.ibm.icu.util.StringTokenizer;

public class TimeSorter implements Comparator<String> {
    public int compare( String a, String b ) {
        int aValue = getUnitValue( a ) + getValue( a );
        int bValue = getUnitValue( b ) + getValue( b );
        return aValue - bValue;
    }

    public int getUnitValue( String a ) {
        if ( a.contains( "yrs" ) ) return 20000;
        if ( a.contains( "mos" ) || a.contains( "mo" )) return 10000;
        if ( a.contains( "pcw" ) ) return 0;
        return 0;
    }

    public int getValue( String a ) {
        StringTokenizer tok = new StringTokenizer( a );
        return Integer.parseInt( tok.nextToken() );
    }

}
