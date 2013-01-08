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
package ubic.BAMSandAllen.geneFilters;

import java.util.List;

import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class ABA1000ConflictsGeneFilter extends PrefixGeneFilter implements GeneFilter {

    public String getName() {
        return "Allen 1000 validation set gene filter";
    }

    public ABA1000ConflictsGeneFilter() {
        super();
        prefixes.add( "Nrp2[" );
        prefixes.add( "Trp53bp1[" );
        prefixes.add( "Ntrk2[" );
        prefixes.add( "Titf1[" );
        prefixes.add( "C920006C10Rik[" );
        prefixes.add( "Cdk4[" );
        prefixes.add( "Wnt7b[" );
        prefixes.add( "Ppfibp1[" );
        prefixes.add( "Galr1[" );
        prefixes.add( "Lhx2[" );
        prefixes.add( "Gad2[" );
        prefixes.add( "Gabrd[" );
        prefixes.add( "Hes1[" );
        prefixes.add( "Cxcl12[" );
        prefixes.add( "Sepm[" );
        prefixes.add( "Unc5d[" );
        prefixes.add( "Tlx3[" );
        prefixes.add( "Mlf1[" );
        prefixes.add( "Fxr1h[" );
        prefixes.add( "Axot[" );
        prefixes.add( "Fgf15[" );
        prefixes.add( "Ceacam10[" );
        prefixes.add( "Dedd[" );
        prefixes.add( "P2rx7[" );
        prefixes.add( "P2ry1[" );
        prefixes.add( "Ranbp9[" );
        prefixes.add( "Wwc1[" );
        prefixes.add( "Crhr2[" );
    }
}
