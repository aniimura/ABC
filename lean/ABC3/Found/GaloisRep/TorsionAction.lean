import ABC3.Found.GaloisRep.WeilBilinFull

/-!
# Galois (G5) 第 196 ブロック —— **★★★★★★★`E[n]` の関数体への作用と Artin**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★非退化性が要求するもの

非退化性(5 性質の最後)の古典的な証明は次の形をしている:

    e_n(P, ·) ≡ 1  ⟹  τ_Q(g) = g (∀Q ∈ E[n])  ⟹  g ∈ F(E)^{E[n]} = [n]^* F(E)
                   ⟹  g = h ∘ [n]  ⟹  h^n = f_P  ⟹  div(h) = (P) − (O)  ⟹  P = O

★最後の段は **`Point.toClass` の単射性**(mathlib)でそのまま済む。
★★止まるのは **`F(E)^{E[n]} = [n]^* F(E)`**——`[n]` の Galois 理論である。

## ★★★★★★★本ブロック——包含と Artin

| 段 | 状態 |
|---|---|
| `E[n]` が `F(E)` に**忠実に**作用する | ✅ 本ブロック |
| `[F(E) : F(E)^{E[n]}] = n²`(Artin) | ✅ 本ブロック |
| `[n]^* F(E) ⊆ F(E)^{E[n]}` | ✅ 本ブロック |
| **`[F(E) : [n]^* F(E)] = n²`**(`deg[n] = n²`) | ❌ **残る** |

★★★等号 `F(E)^{E[n]} = [n]^* F(E)` を出すには最後の 1 行が要る。
★★これは `deg[n] = n²` であり、**分割多項式(`Φ_n/Ψ_n²` の次数が `n²`)**か
**双対同種**のどちらかが要る——本プロジェクトがこれまで避けてきた量である。

### ★忠実性の中身

`translateAut` は `IsTranslate` で一意に決まる(第 190 `isTranslate_unique`)ので、
`τ_Q = τ_{Q'}` なら第 189 の `isTranslate_iff` から `toFF Q = toFF Q'`、
`Point.map_injective` で `Q = Q'`。

### ★Artin の適用

mathlib の `FixedPoints.finrank_eq_card` に `Fintype` と `FaithfulSMul` を与える。
★位数は `#E[n] = n²`(`torsion_card`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `translateAut_mul` / `_one` | `Q ↦ τ_Q` は準同型 |
| `translateAut_injective` | ★★★★★**忠実** |
| `torsionAutHom` | ★★★★★★**`E[n] →* Aut(F(E))`** |
| `card_torsGroup` | `#E[n] = n²` |
| `finrank_fixedPoints_torsGroup` | ★★★★★★★**Artin: `[F(E) : F(E)^{E[n]}] = n²`** |
| `muExt_mem_fixedPoints` | ★★★★★★**`[n]^* F(E) ⊆ F(E)^{E[n]}`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine IsDedekindDomain

variable {F : Type} [Field F] [DecidableEq F] [IsAlgClosed F] [Infinite F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic] [inst : IsDedekindDomain W.CoordinateRing]

/-! ## ★★★★★平行移動は群準同型で忠実 -/

theorem translateAut_mul (h4 : (4 : F) ≠ 0) (Q₁ Q₂ : W.Point) :
    translateAut W h4 Q₁ * translateAut W h4 Q₂ = translateAut W h4 (Q₁ + Q₂) := by
  rw [← translateAut_trans W h4 Q₁ Q₂]
  rfl

theorem translateAut_one (h4 : (4 : F) ≠ 0) : translateAut W h4 (0 : W.Point) = 1 :=
  isTranslate_translateAut W h4 0

/-- ★★★★★**忠実**——`τ_Q` は `Q` を決める。 -/
theorem translateAut_injective (h4 : (4 : F) ≠ 0) :
    Function.Injective (translateAut W h4) := by
  intro Q Q' hQ
  have h1 := (isTranslate_iff W (translateAut W h4 Q) Q).1 (isTranslate_translateAut W h4 Q)
  have h2 := (isTranslate_iff W (translateAut W h4 Q') Q').1 (isTranslate_translateAut W h4 Q')
  rw [hQ, h2] at h1
  refine WeierstrassCurve.Affine.Point.map_injective (W' := W)
    (Algebra.ofId F W.FunctionField) ?_
  calc toFF W Q
      = -(ABC3.Found.GaloisRep.genericPoint W)
        + (ABC3.Found.GaloisRep.genericPoint W + toFF W Q) := by abel
    _ = -(ABC3.Found.GaloisRep.genericPoint W)
        + (ABC3.Found.GaloisRep.genericPoint W + toFF W Q') := by rw [h1]
    _ = toFF W Q' := by abel

/-! ## ★★★★★★`E[n] →* Aut(F(E))` -/

/-- `E[n]` を乗法的に見た群。 -/
abbrev torsGroup (n : ℕ) := Multiplicative ((nsmulEndo (A := W.Point) n).ker)

/-- ★★★★★★**`E[n]` から関数体の自己同型群への単射準同型**。 -/
noncomputable def torsionAutHom (h4 : (4 : F) ≠ 0) (n : ℕ) :
    torsGroup W n →* (W.FunctionField ≃ₐ[F] W.FunctionField) where
  toFun Q := translateAut W h4 ((Multiplicative.toAdd Q : (nsmulEndo (A := W.Point) n).ker) :
    W.Point)
  map_one' := translateAut_one W h4
  map_mul' _ _ := (translateAut_mul W h4 _ _).symm

theorem torsionAutHom_injective (h4 : (4 : F) ≠ 0) (n : ℕ) :
    Function.Injective (torsionAutHom W h4 n) := by
  intro Q Q' hQ
  exact Multiplicative.toAdd.injective (Subtype.ext (translateAut_injective W h4 hQ))

/-! ## ★位数 -/

theorem finite_torsGroup (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0) : Finite (torsGroup W n) := by
  have hcard := torsion_card W hΔ n hn hchar
  have hfin : Finite {P : W.Point // n • P = 0} :=
    Nat.finite_of_card_ne_zero (by rw [hcard]; positivity)
  have hequiv : ((nsmulEndo (A := W.Point) n).ker) ≃ {P : W.Point // n • P = 0} :=
    Equiv.subtypeEquivRight (fun P => mem_ker_nsmulEndo n P)
  haveI : Finite ((nsmulEndo (A := W.Point) n).ker) := Finite.of_equiv _ hequiv.symm
  exact inferInstanceAs (Finite (Multiplicative _))

theorem card_torsGroup (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0) :
    Nat.card (torsGroup W n) = n ^ 2 := by
  have hequiv : ((nsmulEndo (A := W.Point) n).ker) ≃ {P : W.Point // n • P = 0} :=
    Equiv.subtypeEquivRight (fun P => mem_ker_nsmulEndo n P)
  rw [show Nat.card (torsGroup W n) = Nat.card ((nsmulEndo (A := W.Point) n).ker) from rfl,
    Nat.card_congr hequiv]
  exact torsion_card W hΔ n hn hchar

/-! ## ★★★★★★★Artin -/

/-- ★★★★★★★**Artin: `[F(E) : F(E)^{E[n]}] = n²`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `FixedPoints.finrank_eq_card` に忠実性と位数 `n²` を与えるだけ。 -/
theorem finrank_fixedPoints_torsGroup (h4 : (4 : F) ≠ 0) (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0) :
    letI := MulSemiringAction.compHom W.FunctionField (torsionAutHom W h4 n)
    Module.finrank (FixedPoints.subfield (torsGroup W n) W.FunctionField) W.FunctionField
      = n ^ 2 := by
  haveI := finite_torsGroup W hΔ n hn hchar
  haveI : Fintype (torsGroup W n) := Fintype.ofFinite _
  letI := MulSemiringAction.compHom W.FunctionField (torsionAutHom W h4 n)
  haveI : FaithfulSMul (torsGroup W n) W.FunctionField := by
    refine ⟨fun {g₁ g₂} hsm => ?_⟩
    refine torsionAutHom_injective W h4 n ?_
    exact AlgEquiv.ext (fun x => hsm x)
  rw [FixedPoints.finrank_eq_card (torsGroup W n) W.FunctionField, ← Nat.card_eq_fintype_card]
  exact card_torsGroup W hΔ n hn hchar

/-- ★★★★★★**`[n]^* F(E)` は固定体に入る**。 -/
theorem muExt_mem_fixedPoints (h4 : (4 : F) ≠ 0) (n : ℕ)
    {μ : W.CoordinateRing →+* W.FunctionField} (hinj : Function.Injective μ)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (z : W.FunctionField) :
    letI := MulSemiringAction.compHom W.FunctionField (torsionAutHom W h4 n)
    muExt W hinj z ∈ FixedPoints.subfield (torsGroup W n) W.FunctionField := by
  letI := MulSemiringAction.compHom W.FunctionField (torsionAutHom W h4 n)
  intro Q
  show translateAut W h4 _ (muExt W hinj z) = muExt W hinj z
  refine aut_comp_muExt W (translateAut W h4 _) hinj ?_ z
  refine RingHom.ext (fun r => ?_)
  exact translateAut_fix_mu W h4 n μ hμF hns hμx hμy hμP
    ((mem_ker_nsmulEndo n _).1 (Multiplicative.toAdd Q).2) r

/-! ## ★出典の紐付け(`.src`) -/

def torsionAutHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——E[n] の関数体への作用)",
    sectionId := "genell-thm-3-8" }

def finrank_fixedPoints_torsGroup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——Artin の定理)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
