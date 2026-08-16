import ABC3.Found.FrdI.Def23

/-!
# [FrdI] Proposition 2.5 —— The Unit-linear Frobenius Functor

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.48–p.50。

原文 (FrdI p.48):
> (The Unit-linear Frobenius Functor) Let

## ★この命題の規模(測定)

**3 条、主張は 4**:

| 条 | # | 内容 |
|---|---|---|
| (i) | 1 | `𝒪^▷(A)^char → Φ(A)` が**全単射**(`Proposition 2.2, (iii)` の包含が実は全単射) |
| (ii) | 2 | `𝒞^istr` が **base-trivial 型** |
| (ii) | 3 | `𝒞^istr` のどの対象も **Frobenius-trivial** |
| (iii) | 4 | **unit-linear Frobenius 函手** `Ψ : 𝒞 ≃ 𝒞(d)`((a) 対象と等長射の上で恒等、(b) `d` の Frobenius 函手と 1-compatible) |

★**(iii) は `Definition 2.4` の `Λ` / `d ∈ Λ>0` / `𝒞(d)` を要する**ので、
ここでは (i)(ii) を実装する。★**(iii) が入るまで `.src` は付けない。**

## ★仮定について

原文は `𝒞` が **Frobenius-normalized・metrically trivial・Aut-ample 型**であるとし、
`τ` を characteristic splitting、`Λ` を `Φ` を supports する monoid type、`d ∈ Λ>0` とする。
★**(i) が使うのは metrically trivial と Aut-ample だけ**(原文の証明どおり)、
★**(ii) が使うのは metrically trivial だけ**である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★(i) —— `𝒪^▷(A)^char → Φ(A)` は全単射

原文 (FrdI p.48):
> sition 2.2, (iii), is, in fact, a bijection.

原文 (FrdI p.49):
> metrically trivial and Aut-ample type. Next, we consider assertion (ii). Since C

★**原文の証明どおりの 3 手**:
1. `Definition 1.3, (iii), (d)` の圏同値で、`c ∈ Φ(A)` を
   co-angular pre-step `ψ : A ⟶ X`(`Div ψ = c`)として実現する(`coaPre_realize`)
2. ★**metrically trivial 型**なので `X ≅ A`、すなわち `ψ` は `A` の自己射に化ける
3. ★**Aut-ample 型**なので、ずれた底の自己同型を `𝒞` の自己同型で打ち消せる

★**単射性の側は `Proposition 2.2, (iii)`(`otri_div_eq_iff`)で既に取れている。**
-/

/-- ★★★**[FrdI] Proposition 2.5, (i)** —— `Div : 𝒪^▷(A) → Φ(A)` は**全射**。

`Proposition 2.2, (iii)` の単射性(`otri_div_eq_iff`)と合わせて、
`𝒪^▷(A)^char → Φ(A)` が**全単射**であることを与える。 -/
theorem prop_2_5_i_surjective (G : Frobenioid P) {A : C}
    (hmt : IsMetricallyTrivial P A) (haa : IsAutAmple P A)
    (c : Φ.val (P.toElem.obj A).base) :
    ∃ x : OTri P A, P.Div ((x : End A) : A ⟶ A) = c := by
  -- 手 1: `c` を co-angular pre-step として実現する
  obtain ⟨X, ψ, hψc, hψs, hψd⟩ := coaPre_realize P G A c
  -- 手 2: metrically trivial 型で `X ≅ A`
  obtain ⟨e⟩ := hmt X ψ hψc hψs
  -- 手 3: Aut-ample 型でずれた底を打ち消す
  have hgi : IsIso (P.Base (ψ ≫ e.hom)) := by
    haveI : IsIso (P.Base ψ) := hψs.2
    haveI : IsIso (P.Base e.hom) := isBaseIsomorphism_of_isIso P e.hom
    rw [P.Base_comp]
    infer_instance
  obtain ⟨φ, hφi, hφb⟩ := haa (P.Base (ψ ≫ e.hom)) hgi
  haveI := hφi
  have hbinv : P.Base (ψ ≫ e.hom) ≫ P.Base (inv φ) = 𝟙 _ := by
    rw [← hφb, ← P.Base_comp, IsIso.hom_inv_id, P.Base_id]
  refine ⟨⟨ψ ≫ e.hom ≫ inv φ, ?_, ?_⟩, ?_⟩
  · -- base-identity
    show P.Base (ψ ≫ e.hom ≫ inv φ) = P.Base (𝟙 A)
    rw [P.Base_id, ← Category.assoc, P.Base_comp, hbinv]
  · -- linear
    show P.degFr (ψ ≫ e.hom ≫ inv φ) = 1
    rw [P.degFr_comp, P.degFr_comp, show P.degFr ψ = 1 from hψs.1,
      degFr_of_isIso P e.hom, degFr_of_isIso P (inv φ)]
    simp
  · -- `Div = c`
    show P.Div (ψ ≫ e.hom ≫ inv φ) = c
    haveI : IsIso (e.hom ≫ inv φ) := inferInstance
    rw [P.Div_comp, show P.Div (e.hom ≫ inv φ) = 0 from isIsometric_of_isIso P _,
      degFr_of_isIso P (e.hom ≫ inv φ), map_zero]
    simpa using hψd

/-! ## ★(ii) —— `𝒞^istr` は base-trivial 型で、どの対象も Frobenius-trivial

原文 (FrdI p.48):
> (ii) Cistr is of base-trivial type. Moreover, every object of Cistr is Frobenius-

原文 (FrdI p.49):
> follows from the existence of Frobenius-trivial objects [cf. Definition 1.3, (i), (a)]

★**原文の証明どおり**:
1. `Definition 1.3, (i), (b)`(`preStepSpan`)で、底が同型な `𝒞^istr` の 2 対象は
   pre-step の span で結ばれる。★**span の頂点も `𝒞^istr` の対象**なので
   `Proposition 1.4, (i)` により両辺の pre-step は自動で co-angular
2. ★**metrically trivial 型**なので、その span の両端は頂点と同型 —— よって互いに同型
3. `Definition 1.3, (i), (a)` を **`𝒞^istr` 自身**に当てる
   (`Proposition 1.9, (v)`、`istr_frobenioidCore`)と Frobenius-trivial 対象が取れ、
   2 の base-trivial 性で**すべての対象へ移る**
-/

variable (F : FrobenioidCore P)

/-- ★★★**[FrdI] Proposition 2.5, (ii) の前半** —— `𝒞^istr` は **base-trivial 型**。 -/
theorem prop_2_5_ii_baseTrivial (hmt : ∀ A : C, IsMetricallyTrivial P A) (A : Istr P) :
    IsBaseTrivial (istrPre P F) A := by
  intro Dd hbi
  obtain ⟨α⟩ := hbi
  obtain ⟨X, σ, τ, hσ, hτ, -⟩ :=
    (istr_frobenioidCore P F).preStepSpan A Dd α.hom inferInstance
  obtain ⟨eA⟩ := hmt X.obj A.obj σ.hom
    (isCoAngular_of_isotropic_dom P F X.property σ.hom) hσ
  obtain ⟨eD⟩ := hmt X.obj Dd.obj τ.hom
    (isCoAngular_of_isotropic_dom P F X.property τ.hom) hτ
  exact ⟨ObjectProperty.isoMk _ (eD ≪≫ eA.symm)⟩

/-- ★★★**[FrdI] Proposition 2.5, (ii) の後半** —— `𝒞^istr` のどの対象も
**Frobenius-trivial**。

★`Definition 1.3, (i), (a)` を `𝒞^istr` 自身に当て、前半の base-trivial 性で移す。 -/
theorem prop_2_5_ii_frobTrivial (hmt : ∀ A : C, IsMetricallyTrivial P A) (A : Istr P) :
    IsFrobeniusTrivial (istrPre P F) A := by
  obtain ⟨B, hB, ⟨e⟩⟩ := (istr_frobenioidCore P F).baseSurj
    ((istrPre P F).toElem.obj A).base
  obtain ⟨θ⟩ := prop_2_5_ii_baseTrivial P F hmt A B ⟨e.symm⟩
  exact isFrobeniusTrivial_of_iso (istrPre P F) (istr_frobenioidCore P F) θ hB

end ABC3.Found.FrdI
