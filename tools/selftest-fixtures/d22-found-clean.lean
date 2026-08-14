-- D22 fixture: Found/ が sorry を含まない(通るべき)
-- docstring 内で "sorry" に言及しても誤検出しないことも兼ねて確認する。
namespace Fixture

/-- この定理は sorry を使っていない。 -/
theorem complete : (0 : Nat) = 0 := rfl

end Fixture
