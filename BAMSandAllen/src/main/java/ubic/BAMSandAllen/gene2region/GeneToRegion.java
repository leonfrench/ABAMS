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
package ubic.BAMSandAllen.gene2region;


/**
 * Simple class to store gene to region correlations
 * 
 * @author leon
 *
 */
public class GeneToRegion implements Comparable<GeneToRegion> {
    public String gene;
    public String region;
    public double correlation;

    public GeneToRegion( String gene, String region, double correlation ) {
        super();
        this.gene = gene;
        this.region = region;
        this.correlation = correlation;
    }

    public int compareTo( GeneToRegion g2r ) {
        return Double.compare( this.correlation, g2r.correlation );
    }

    public String toString() {
        return gene + " vrs " + region + " cor:" + correlation;
    }
}
