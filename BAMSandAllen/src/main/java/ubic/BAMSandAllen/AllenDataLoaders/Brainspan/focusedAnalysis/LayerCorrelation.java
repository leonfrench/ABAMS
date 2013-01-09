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
package ubic.BAMSandAllen.AllenDataLoaders.Brainspan.focusedAnalysis;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import cern.colt.list.DoubleArrayList;

import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.Brainspan.BrainSpanGeneric;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.dataStructure.params.ParamKeeper;
import ubic.basecode.math.CorrelationStats;
import ubic.basecode.math.DescriptiveWithMissing;

public class LayerCorrelation extends BrainSpanGeneric {
    private static Log log = LogFactory.getLog( BrainSpanGeneric.class.getName() );
    String age;
    Pattern pattern;

    public LayerCorrelation( String age ) throws Exception {
        super( "blueprint" );
        this.age = age;

        String patternStr = "ayer.([0-6])";

        // Compile and use regular expression
        pattern = Pattern.compile( patternStr );

        Set<String> regions = getUniqueRegions();
        // DoubleMatrix<String, String> layerData = new DenseDoubleMatrix<String, String>( uniqueGenes.size(), 5 );
        // DoubleMatrix<String, String> layerData = new DenseDoubleMatrix<String, String>( 1, 5 );
        List<String> layerRegions = new LinkedList<String>();
        List<Double> layerValues = new LinkedList<Double>();

        for ( String region : regions ) {
            String layer = getLayer( region );
            if ( layer != null ) {
                log.info( region + " -> " + layer );
                // layerData.addColumnName( region );
                layerRegions.add( region );
                layerValues.add( Double.parseDouble( layer ) );
            }
        }

        double[] layerValueArray = ArrayUtils.toPrimitive( layerValues.toArray( new Double[0] ) );

        Set<String> uniqueGenes = getUniqueGenes();
        // layerData.addColumnName( "2" );
        // layerData.addColumnName( "3" );
        // layerData.addColumnName( "4" );
        // layerData.addColumnName( "5" );
        // layerData.addColumnName( "6" );

        ParamKeeper keeper = new ParamKeeper();

        int geneCount = 0;
        for ( String gene : uniqueGenes ) {
            if ( ++geneCount % 10 == 0 ) log.info( "Gene count:" + geneCount );
            // if ( geneCount > 200 ) break;

            HashMap<String, String> line = new HashMap<String, String>();
            line.put( "gene", gene );

            List<Double> geneValues = new LinkedList<Double>();

            Map<String, List<Double>> dataMap = findData( gene, age, layerRegions );
            for ( String region : layerRegions ) {
                List<Double> data = dataMap.get( region );
                double[] valueArray = ArrayUtils.toPrimitive( data.toArray( new Double[0] ) );
                DoubleArrayList expValuesDAL = new DoubleArrayList( valueArray );
                double mean = DescriptiveWithMissing.mean( expValuesDAL );
                geneValues.add( mean );
            }

            double[] geneValueArray = ArrayUtils.toPrimitive( geneValues.toArray( new Double[0] ) );
            DoubleArrayList expValuesDAL = new DoubleArrayList( geneValueArray );
            double averageExp = DescriptiveWithMissing.mean( expValuesDAL );

            double spearmanCor = Util.spearmanCorrel( geneValueArray, layerValueArray );
            double spearmanP = CorrelationStats.spearmanPvalue( spearmanCor, geneValueArray.length )
                    * uniqueGenes.size();
            double pearsonCor = CorrelationStats.correl( geneValueArray, layerValueArray );
            double pearsonP = CorrelationStats.pvalue( pearsonCor, layerValueArray.length ) * uniqueGenes.size();

            line.put( "variance", "" + Util.variance( geneValueArray ) );
            line.put( "spearmanCor", "" + spearmanCor );
            line.put( "spearmanP", "" + spearmanP );
            line.put( "pearsonCor", "" + pearsonCor );
            line.put( "pearsonP", "" + pearsonP );
            line.put( "averageExp", "" + averageExp );
            keeper.addParamInstance( line );
        }
        keeper.writeExcel( "/tmp/layercor." + age + ".xls" );
    }

    public String getLayer( String name ) {
        Matcher matcher = pattern.matcher( name );

        if ( matcher.find() ) {
            return matcher.group( 1 );
        }
        return null;

    }

    public static void main( String[] args ) throws Exception {
        LayerCorrelation loader;
        loader = new LayerCorrelation( "3 mo" );
        loader = new LayerCorrelation( "12 mo" );
        loader = new LayerCorrelation( "48 mo" );
    }

}
