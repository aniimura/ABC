import ABC3.Skeleton.PGC.Setup

/-!
# [pGC] Corollary 1.3 — 素朴な Interface 設計が退化することの実証

Corollary 1.3 は「惰性群 `I_K ⊆ Γ_K` が Γ_K から群論的に決定できる」と主張する。

素朴には、`I_K` を **`Interface` の自由なデータ**として置きたくなる:

```lean
structure InertiaData where
  inertia : (K : PAdicLocalField p) → Subgroup K.absGal   -- 条件なし
```

このファイルは、**その設計では Corollary 1.3 が自明に成り立ってしまう**ことを
機械で示す。`I_K := ⊥` と置けば通る。

## なぜこれが重要か

これは「偽になる」のではなく「**自明になる**」型の失敗である。
`lake build` も `sorry` 検査も `#print axioms` も通る。
`ABC3/Meta/Calibration.lean` の実演と同型だが、そこでは *仮説* が空虚だったのに対し、
ここでは *結論* が退化した witness で満たされる——**別の顔をした同じ病**。

## 結論

`I_K` は **posit するものではなく derive するもの**である。
原文 p.3 の不分岐判定(「L is unramified over K if and only q_L = q^[Γ_K:H]」)から
構成しなければ、Corollary 1.3 は内容を持たない。
その構成は原文が明示していない(暗黙の段 #7)。

**これは原典の主張ではない**(我々のモデルについての事実)ので `.src` を持たない。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Meta

variable {p : ℕ} [Fact p.Prime]

/-- 惰性群を「条件なしの自由なデータ」として置いた場合の、K に付随する対象。

移送は部分群の押し出し `S ↦ α(S)`——原典の「α によって対応物へ移る」を
部分群について読めばこれになる。 -/
noncomputable def inertiaObjectOf (I : (K : PAdicLocalField p) → Subgroup K.absGal) :
    AssociatedObject p where
  Obj := fun K => Subgroup K.absGal
  obj := I
  transport := fun α S => S.map α.toMulEquiv.toMonoidHom

/-- ★**退化の実証**: `I_K := ⊥` と置くと Corollary 1.3 の形は成り立ってしまう。

自明群は任意の群準同型で自明群へ移るので、「どの同型で写しても対応が保たれる」が
無内容に成り立つ。`sorry` 無し。 -/
theorem bot_inertia_recoverable :
    (inertiaObjectOf (p := p) (fun _ => ⊥)).RecoverableFromAbsGal := by
  intro _ _ _
  exact Subgroup.map_bot _

/-- 同様に `I_K := ⊤` でも成り立つ——退化は片側だけではない。 -/
theorem top_inertia_recoverable :
    (inertiaObjectOf (p := p) (fun _ => ⊤)).RecoverableFromAbsGal := by
  intro _ _ α
  exact Subgroup.map_top_of_surjective _ α.toMulEquiv.surjective

-- 受理ゲートを全部通ることの確認(退化した witness であるにもかかわらず)
#print axioms bot_inertia_recoverable
#print axioms top_inertia_recoverable

/-!
## 読み

`⊥` と `⊤` という**両極の退化 witness**が、どちらも Corollary 1.3 の形を満たす。
ゆえに「`I_K` を自由なデータとして置き、`RecoverableFromAbsGal` を要求する」という
設計は、惰性群について何も言っていない。

正しい設計は、`I_K` を原文の不分岐判定から**構成する**こと:

```
I_K = ⋂ { H ⊆ Γ_K : H は開、かつ L_H/K が不分岐 }
```

そして原文の判定条件は `q_{L_H} = q^{[Γ_K:H]}`。これを書くには、
「開部分群 H に対応する中間体 L_H もまた p進局所体であり、その絶対 Galois 群が H」
という対応が要る——原文 p.3 が Proposition 1.2 を (L, H) に適用するときに
前提しているもの。ここが次の `Interface` になる。
-/

end ABC3.Check.PGC
