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

import ubic.BAMSandAllen.SetupParameters;

public class NewStructureExpressionInformationLoader extends AllenCSVLoader {
    public static final int imageIDPosition = 0;
    public static final int densityPosition = 2;
    public static final int levelPosition = 3;
    public static final int energyPosition = 4;
    ImageSeriesInfoLoader infoLoader;

    public NewStructureExpressionInformationLoader() throws Exception {
        super();
        headerLine = true;
        sep = ',';
        filename = SetupParameters.config.getString( "abams.allen.newExpressionData" );
        regionNamePosition = 1;

        infoLoader = new ImageSeriesInfoLoader();
        init();
    }

    public String getRowName( String line[] ) {
        String imageid = line[imageIDPosition];
        String geneShort = infoLoader.getGeneFromImageID( imageid );
        return geneShort + "[" + imageid + "]";
    }

    public static void main( String[] args ) throws Exception {
        NewStructureExpressionInformationLoader loader = new NewStructureExpressionInformationLoader();
        AllenMinerLoader miner = new AllenMinerLoader();
    }

}
