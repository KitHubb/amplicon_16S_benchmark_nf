class TrimmingConfig {
    static Map resolve(params) {
        def result = [:]
        ['trim_poly_a', 'trim_poly_g', 'discard_untrimmed'].each { key ->
            if (!(params[key].toString() in ['true', 'false'])) {
                throw new IllegalArgumentException("${key} must be true or false")
            }
            result[key] = params[key].toString().toBoolean()
        }
        ['extra_adapters_f', 'extra_adapters_r'].each { key ->
            if (!(params[key] instanceof List)) {
                throw new IllegalArgumentException("${key} must be a YAML list of DNA adapter sequences")
            }
            result[key] = params[key].collect { sequence ->
                if (!(sequence instanceof CharSequence) || !(sequence ==~ /[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+/)) {
                    throw new IllegalArgumentException("${key} must contain non-empty IUPAC DNA sequences")
                }
                sequence.toString().toUpperCase()
            }.unique()
        }
        ['poly_min_length', 'min_length', 'quality'].each { key ->
            if (!(params[key].toString() ==~ /\d+/) || params[key].toString().toLong() < (key == 'quality' ? 0 : 1) || params[key].toString().toLong() > Integer.MAX_VALUE) {
                throw new IllegalArgumentException("${key} must be an integer >= ${key == 'quality' ? 0 : 1}")
            }
            result[key] = params[key].toString().toInteger()
        }
        return result
    }
}
