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

import java.io.Serializable;

public class LOOCorrelObject implements Serializable {
    double sxy = 0.0, sxx = 0.0, syy = 0.0, sx = 0.0, sy = 0.0;
    int numused = 0;

    public LOOCorrelObject clone() {
        LOOCorrelObject result = new LOOCorrelObject();
        result.sx = sx;
        result.sy = sy;
        result.sxy = sxy;
        result.sxx = sxx;
        result.numused = numused;
        return result;

    }

    public double correlAll( double[] ival, double[] jval ) {
        int length = Math.min( ival.length, jval.length );
        for ( int k = 0; k < length; k++ ) {
            double xj = ival[k];
            double yj = jval[k];
            addPair( xj, yj );
        }
        return correl();
    }

    public double correl() {
        /* do it the old fashioned way */
        if ( numused < 2 ) return Double.NaN;
        double denom = ( sxx - sx * sx / numused ) * ( syy - sy * sy / numused );
        if ( denom <= 0 ) return Double.NaN;
        double correl = ( sxy - sx * sy / numused ) / Math.sqrt( denom );
        return correl;
    }

    public void removePair( double xj, double yj ) {
        if ( !Double.isNaN( xj ) && !Double.isNaN( yj ) ) {
            // none to remove
            if ( numused == 0 ) throw new RuntimeException();
            sx -= xj;
            sy -= yj;
            sxy -= xj * yj;
            sxx -= xj * xj;
            syy -= yj * yj;
            numused--;
        }
    }

    public void addPair( double xj, double yj ) {
        if ( !Double.isNaN( xj ) && !Double.isNaN( yj ) ) {
            sx += xj;
            sy += yj;
            sxy += xj * yj;
            sxx += xj * xj;
            syy += yj * yj;
            numused++;
        }
    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        // TODO Auto-generated method stub

    }

}
