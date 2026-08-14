import ABC3.Skeleton.PGC.Section1

/-!
# [pGC] §1 — `RecoverableFromAbsGal` の識別力の検査

`ABC3/Skeleton/PGC/Section1.lean` で切り出した「group-theoretically に回復できる」の
定義について、**述語が常に真でも常に偽でもない**ことを確かめる。

なぜ要るか: 原典が暗黙にしていた「移送の規則」(`AssociatedObject.transport`)を
我々がデータとして明示した以上、その選び方次第で述語が空虚に潰れうる——
`ABC3/Meta/Calibration.lean` が実演した失敗と同型。

**これは原典の主張ではない**(我々のモデルについての事実)ので `.src` を持たない。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Meta

variable {p : ℕ} [Fact p.Prime]

/-- `PAdicLocalField p` は空でない(`ℚ_[p]` 自身)。以下の検査の土台。 -/
noncomputable def base : PAdicLocalField p := { carrier := ℚ_[p] }

/-! ## 常に真ではない -/

/-- 移送を反転させた対象。恒等射を取れば `!true = true` を要求することになる。 -/
def flipped : AssociatedObject p where
  Obj := fun _ => Bool
  obj := fun _ => true
  transport := fun _ b => !b

theorem flipped_not_recoverable : ¬ (flipped (p := p)).RecoverableFromAbsGal := by
  intro h
  exact Bool.noConfusion (h (K := base) (K' := base) (ContinuousMulEquiv.refl _))

/-! ## 常に偽でもない -/

/-- 自明な対象は回復可能。 -/
def trivialObj : AssociatedObject p where
  Obj := fun _ => Unit
  obj := fun _ => ()
  transport := fun _ u => u

theorem trivialObj_recoverable : (trivialObj (p := p)).RecoverableFromAbsGal := by
  intro _ _ _; rfl

end ABC3.Check.PGC
