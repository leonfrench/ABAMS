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

import java.util.HashSet;
import java.util.Set;

public class NomenClatureEntry implements Cloneable {

    public NomenClatureEntry() {
        expressedGenes = new HashSet<String>();
    }

    public String acro;
    public String name;
    public String NNID;
    public String source;
    public String NNName;
    public Set<String> expressedGenes;

    public NomenClatureEntry copy() {
        try {
            return ( NomenClatureEntry ) this.clone();
        } catch ( Exception e ) {
            e.printStackTrace();
            return null;
        }
    }
}
