/*
 * MAF indexer based on the normal indexer
 * The index array is sorted by start > end order
 * When a query is sent, instead of the closest upstream key a backwards search
 * is performed starting at the right end of the query slice and if there is an overlap
 * all hits are returned. range up to a limit is searched.
 * very unlikely large overlaps are ignored
 * instead of a fixed step size, a linear search is performed to find the starting index.
 */
