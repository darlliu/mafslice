/*
 * A base indexer for sequentially indexed contents
 * contains a map from string (by chromosome for example) to vector of indices
 * and a sequentially sorted vector with a fixed step size from indices to keys
 * and a on disk hash table from keys to actual content (chunks of sequences)
 * Functionalities are provided to query sequences from any slice of indices on major key
 */


