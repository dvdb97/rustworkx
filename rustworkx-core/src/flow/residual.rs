

#[derive(Clone, Copy)]
pub enum ResArcType<E> {
    Forward(E),
    Backward(E),
}
