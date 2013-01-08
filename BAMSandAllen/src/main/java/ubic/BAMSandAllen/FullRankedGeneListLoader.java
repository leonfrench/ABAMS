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

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.basecode.util.FileTools;

/**
 * A ranked genelistloader that assumes we are using all the genes as opposed to a subset.
 * 
 * @author leon
 */
public class FullRankedGeneListLoader extends RankedGeneListLoader {
    private static Log log = LogFactory.getLog( FullRankedGeneListLoader.class.getName() );

    boolean hasNonExp;

    public FullRankedGeneListLoader( String filename, boolean hasNonExp ) throws Exception {
        super( filename );

        this.hasNonExp = hasNonExp;

        if ( hasNonExp ) addAllenNonExpToEnd();
        addMissingToTop();
    }

    public void removeAllenNonExp() throws Exception {
        List<String> nonExps = getNonExpRows();
        lines.removeAll( nonExps );
    }

    public void addAllenNonExpToEnd() throws Exception {
        List<String> nonExps = getNonExpRows();

        Collections.shuffle( nonExps, new Random( 1 ) );
        log.info( "Non Expressor size:" + nonExps.size() );
        // add to end of list
        for ( String nonExp : nonExps ) {
            // it should not contain it
            if ( !lines.contains( nonExp ) ) lines.add( nonExp );
        }
    }

    private List<String> getNonExpRows() throws Exception {
        List<String> originalRows = getOriginalRows();
        List<String> nonExps = Util.getRowNamesFromGeneNames( SetupParameters.getDataFolder() + "ABANonexpressed.txt",
                originalRows );
        return nonExps;
    }


    public List<String> addMissingToTop() throws Exception {
        return addMissing( true );
    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        // TODO Auto-generated method stub

    }

}
