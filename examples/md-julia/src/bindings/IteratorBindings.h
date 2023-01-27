struct WrapIteratorInterface {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrapedT = typename T::type;        
    }
};

struct WrapIteratorWrapper {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrappedT = typename T::type;
        wrapped.method("isValid", &WrappedT::isValid);
        wrapped.method("inc", &WrappedT::operator++);
        wrapped.method("deref", &WrappedT::operator*);
    }
};