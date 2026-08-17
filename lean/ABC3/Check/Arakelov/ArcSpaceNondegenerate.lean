import ABC3.Found.Arakelov.ArcSpaceInterface

/-!
# 負の対照 —— **`ArcSpaceData` が退化 witness を弾くこと**(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

## ★★★何を確かめるか

C1(`ArcSpaceData`)を書いたあと、**退化 witness を実際に書いて落ちることを確認する**。
★★2026-08-17 の実測で **1 つ穴が見つかった**:

| 退化 witness | 塞いだ要求 |
|---|---|
| `Arc := PUnit`(1 点) | `equivComplexPoints` |
| `topology := ⊥`(離散) | `topology_affine`(★本ファイルの `arcTopologyAffine_ne_bot`) |
| `topology := ⊤`(密着)+ `evalAffine := 0` | ★★★**`evalAffine_spec`**(後から足した) |

★★★**3 番目は当初通っていた。** `evalAffine := fun _ _ _ => 0` とすると
`induced (fun _ _ => 0) Pi = ⊤` なので `topology_affine` は `topology = ⊤` を
要求するだけになり、`topology_openImmersion` も `induced g ⊤ = ⊤` で自明に成り立つ。
★これを塞ぐために `evalAffine_spec`(評価は `Spec.preimage` が与える環準同型)を足した。

## ★★機構

`ℂ[X]` を取ると、各 `c : ℂ` に対し `aeval c : ℂ[X] → ℂ` が点 `ptC c` を与える。
★`c ↦ ptC c` は**連続**(多項式の値は代入点について連続)で**単射**なので、
`{ptC 0}` が開なら `{0} ⊆ ℂ` が開になってしまう——これは偽である。
-/

universe u v

namespace ABC3.Check.Arakelov

open AlgebraicGeometry CategoryTheory TopologicalSpace
open ABC3.Interface.Arakelov ABC3.Found.Arakelov

/-! ## ★★退化構造(固定する要求を外したもの) -/

/-- ★**`equivComplexPoints` / `evalAffine_spec` / `topology_affine` /
`topology_openImmersion` を外した弱い構造**。

★★これが退化 witness を持つことが、外した 4 要求が**効いている**ことの証拠である。 -/
structure ArcSpaceDataWeak where
  Arc : Scheme.{0} → Type
  topology : (X : Scheme.{0}) → TopologicalSpace (Arc X)
  conj : (X : Scheme.{0}) → Arc X → Arc X
  conj_involutive : ∀ (X : Scheme.{0}) (p : Arc X), conj X (conj X p) = p
  conj_continuous : ∀ (X : Scheme.{0}), @Continuous _ _ (topology X) (topology X) (conj X)
  evalAffine : (A : CommRingCat.{0}) → Arc (Spec A) → A → ℂ
  map : {X Y : Scheme.{0}} → (X ⟶ Y) → Arc X → Arc Y

/-- ★★★**1 点集合は弱い構造を満たす**——これが `equivComplexPoints` の存在理由である。 -/
def trivialArcData : ArcSpaceDataWeak where
  Arc := fun _ => PUnit
  topology := fun _ => ⊥
  conj := fun _ p => p
  conj_involutive := fun _ _ => rfl
  conj_continuous := fun _ => continuous_bot
  evalAffine := fun _ _ _ => 0
  map := fun _ p => p

/-- ★★★**構造だけでは 1 点集合が通る**(負の対照)。 -/
theorem structure_alone_permits_trivial : Nonempty ArcSpaceDataWeak :=
  ⟨trivialArcData⟩

/-! ## ★★★`ℂ[X]` の複素点(離散でないことの種) -/

/-- ★`c : ℂ` が定める `Spec ℂ[X]` の複素点(= `X ↦ c` の代入)。 -/
noncomputable def ptC (c : ℂ) :
    Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of (Polynomial ℂ)) :=
  Spec.map (CommRingCat.ofHom (Polynomial.aeval c).toRingHom)

/-- ★`ptC c` での切断 `a` の値は `a.eval c` である。 -/
@[simp] theorem evalAffine_ptC (c : ℂ) (a : Polynomial ℂ) :
    evalAffine (CommRingCat.of (Polynomial ℂ)) (ptC c) a = a.eval c := by
  show (evalHom (CommRingCat.of (Polynomial ℂ)) (ptC c)).hom a = _
  rw [ptC, evalHom_Spec_map]
  show (Polynomial.aeval c) a = _
  simp

/-- ★`c ↦ ptC c` は単射である(`X` を評価すれば `c` が出る)。 -/
theorem ptC_injective : Function.Injective ptC := by
  intro c c' h
  have hx := congrArg
    (fun p => evalAffine (CommRingCat.of (Polynomial ℂ)) p (Polynomial.X)) h
  simpa using hx

/-- ★★`c ↦ ptC c` は `arcTopologyAffine` について連続である。

★機構は「多項式の値は代入点について連続」(`Polynomial.continuous`)+
各点収束位相(`continuous_pi`)。 -/
theorem continuous_ptC :
    @Continuous _ _ _ (arcTopologyAffine (CommRingCat.of (Polynomial ℂ))) ptC := by
  refine continuous_induced_rng.2 (continuous_pi fun a => ?_)
  show Continuous fun c : ℂ => evalAffine (CommRingCat.of (Polynomial ℂ)) (ptC c) a
  simp only [evalAffine_ptC]
  exact a.continuous

/-! ## ★★★★離散位相は `topology_affine` を満たさない -/

/-- ★`{0}` は `ℂ` の開集合ではない。 -/
theorem not_isOpen_singleton_zero : ¬ IsOpen ({0} : Set ℂ) := by
  intro h
  obtain ⟨ε, hε, hball⟩ := Metric.isOpen_iff.1 h 0 rfl
  have hmem : ((ε / 2 : ℝ) : ℂ) ∈ Metric.ball (0 : ℂ) ε := by
    simp only [Metric.mem_ball, dist_zero_right, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos (by positivity : (0:ℝ) < ε / 2)]
    linarith
  have hz := hball hmem
  simp only [Set.mem_singleton_iff, Complex.ofReal_eq_zero] at hz
  linarith

/-- ★★★★**`X^arc` の位相は離散ではない**(負の対照の核心)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `topology_affine` を課す理由である——課さなければ
`topology := ⊥` で `conj_continuous` が自明に通ってしまう。 -/
theorem arcTopologyAffine_ne_bot :
    (⊥ : TopologicalSpace (Spec (CommRingCat.of ℂ) ⟶
        Spec (CommRingCat.of (Polynomial ℂ))))
      ≠ arcTopologyAffine (CommRingCat.of (Polynomial ℂ)) := by
  letI := arcTopologyAffine (CommRingCat.of (Polynomial ℂ))
  intro h
  have hopen : @IsOpen _ (arcTopologyAffine (CommRingCat.of (Polynomial ℂ)))
      {ptC 0} := by
    rw [← h]; trivial
  refine not_isOpen_singleton_zero ?_
  have hpre : ptC ⁻¹' {ptC 0} = ({0} : Set ℂ) := by
    ext c
    simp only [Set.mem_preimage, Set.mem_singleton_iff]
    exact ⟨fun hc => ptC_injective hc, fun hc => by rw [hc]⟩
  rw [← hpre]
  exact continuous_ptC.isOpen_preimage _ hopen

/-- ★★★**実装の位相は離散ではない**(前定理の帰結)。 -/
theorem arcSpaceDataImpl_topology_ne_bot :
    arcSpaceDataImpl.topology (Spec (CommRingCat.of (Polynomial ℂ))) ≠ ⊥ := fun h =>
  arcTopologyAffine_ne_bot (h.symm.trans (arcTopology_spec _))

/-! ## ★★★★密着位相 + 零評価も落ちる -/

/-- ★★`arcBasicOpen ℂ[X] X` は空でも全体でもない。 -/
theorem arcBasicOpen_X_proper :
    ptC 1 ∈ arcBasicOpen (CommRingCat.of (Polynomial ℂ)) Polynomial.X ∧
      ptC 0 ∉ arcBasicOpen (CommRingCat.of (Polynomial ℂ)) Polynomial.X := by
  constructor
  · show evalAffine _ (ptC 1) Polynomial.X ≠ 0
    simp
  · show ¬ (evalAffine _ (ptC 0) Polynomial.X ≠ 0)
    simp

/-- ★★★★**`X^arc` の位相は密着でもない**。

★★★これが `evalAffine_spec` を課す理由である——課さなければ
`evalAffine := 0` と `topology := ⊤` の組が 7 要求すべてを通ってしまう
(`induced (fun _ _ => 0) Pi = ⊤` かつ `induced g ⊤ = ⊤`)。 -/
theorem arcTopologyAffine_ne_top :
    (⊤ : TopologicalSpace (Spec (CommRingCat.of ℂ) ⟶
        Spec (CommRingCat.of (Polynomial ℂ))))
      ≠ arcTopologyAffine (CommRingCat.of (Polynomial ℂ)) := by
  intro h
  obtain ⟨h1, h0⟩ := arcBasicOpen_X_proper
  have hopen : @IsOpen _ ⊤ (arcBasicOpen (CommRingCat.of (Polynomial ℂ)) Polynomial.X) := by
    rw [h]
    exact isOpen_arcBasicOpen _ _
  rcases (TopologicalSpace.isOpen_top_iff _).1 hopen with he | he
  · rw [he] at h1; exact h1
  · rw [he] at h0; exact h0 (Set.mem_univ _)

/-- ★★★**実装の位相は密着でもない**。 -/
theorem arcSpaceDataImpl_topology_ne_top :
    arcSpaceDataImpl.topology (Spec (CommRingCat.of (Polynomial ℂ))) ≠ ⊤ := fun h =>
  arcTopologyAffine_ne_top (h.symm.trans (arcTopology_spec _))

end ABC3.Check.Arakelov
