/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Ore

/-!
# [FrdI] `Proposition 4.4, (ii)` —— 零因子準同型の**核**

★`ResearchPaper/frdi-decomposition.json` の `otricomm` チェーンの層 1、葉 `ker-eq-image`。

## ★何を示すか

`Prop44Ore.lean` で `𝒪^▷(A^birat) → Φ^gp(A)`(`otriDivGpHom`)を作った。
このファイルはその**核**を決める:

  ★★★**`biratDivGp u = 0` ⟺ `u` は `𝒪^×(A)_𝒞` の像**(`otriDivGp_eq_zero_iff`)

## ★証明の骨(2 行で言うと)

`u = [c]⁻¹ ≫ [f]`(`c` は co-angular pre-step、`Base f = Base c`、`deg f = 1`)と書くと

  `biratDivGp u = (Base c)^{-1,*}(Div f - Div c)`

である。`Φ(A)` は **divisorial**(`PreFrobenioid.divisorial`)なので **integral**、すなわち
`toGp` は単射。ゆえに `biratDivGp u = 0` ⟺ **`Div f = Div c`**。
このとき `c` と `f` は base-equivalent かつ metrically equivalent な co-angular pre-step
なので、`Definition 1.3, (vi)`(`faithfulUpToUnits`)が
`α ∈ 𝒪^×(A)` で `f = c ≫ α` を与える。ゆえに `u = [c]⁻¹ ≫ [c] ≫ [α] = [α]`。

★逆向きは `Φ(A)` が **sharp** であること —— `α` が可逆なら `Div α + Div α⁻¹ = 0` なので
`Div α` は可逆元、よって `0`(`div_eq_zero_of_mem_otimes`)。

★★**divisorial の 2 つの条(integral と sharp)が、両方向でそれぞれ 1 回ずつ効いている。**

## ★これで何が言えるようになるか

★`otricomm` チェーンの層 2 `central-ext` —— **`𝒪^×(A^birat)` は中心拡大**

  `1 → 𝒪^×(A)_𝒞 の像 → 𝒪^×(A^birat) → Φ^gp(A) の部分群 → 1`

であり、核は中心に入る(2026-08-17 に取得済み)。★★したがって可換性の問題は
**交換子が定める交代双線形形式が消えるか**に落ちる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2

section GpInj

variable {D : Type u} [Category.{v} D] {Φ : MonoidOn.{v, u, w} D}

/-- ★`g` が同型なら `gpMap (Φ.map g)` は単射。

★`Φ` は反変なので `Φ.map (inv g)` が逆を与える —— **逆写像を選ぶ必要はない**
(`gpMap_phi_comp` と `MonoidOn.map_id` だけで済む)。 -/
theorem gpMap_map_injective_of_isIso {X Y : D} (g : X ⟶ Y) [IsIso g] :
    Function.Injective (gpMap (Φ.val Y) (Φ.map g)) := by
  intro a b hab
  have key : ∀ z : Gp (Φ.val Y),
      gpMap (Φ.val X) (Φ.map (inv g)) (gpMap (Φ.val Y) (Φ.map g) z) = z := by
    intro z
    rw [gpMap_phi_comp Φ (inv g) g z, IsIso.inv_hom_id]
    have h : (Φ.map (𝟙 Y) : Φ.val Y →+ Φ.val Y) = AddMonoidHom.id _ := by
      ext x; exact MonoidOn.map_id Φ Y x
    rw [h, gpMap_id]
    rfl
  rw [← key a, ← key b, hab]

end GpInj

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {G : Frobenioid P}

/-! ## ★1. `𝒪^×(A)` の元の零因子は `0`(sharp が効く) -/

/-- ★**`𝒪^×(A)` の元は零因子 `0` を持つ**。

`α ≫ α⁻¹ = 𝟙` から `Div α + Div α⁻¹ = 0`、すなわち `Div α` は `Φ(A)` の可逆元。
★`Φ(A)` は **divisorial**、特に **sharp** なので可逆元は `0` だけである。 -/
theorem div_eq_zero_of_mem_otimes {A : C} {α : End A} (hα : α ∈ OTimes P A) :
    P.Div ((α : A ⟶ A)) = 0 := by
  obtain ⟨⟨hb, hl⟩, hu⟩ := hα
  obtain ⟨v, hv⟩ := hu
  subst hv
  have hvv : ((v.val : End A) : A ⟶ A) ≫ ((v.inv : End A) : A ⟶ A) = 𝟙 A := v.inv_val
  have hbase : P.Base ((v.val : End A) : A ⟶ A) = 𝟙 _ := by
    rw [show P.Base ((v.val : End A) : A ⟶ A) = P.Base (𝟙 A) from hb, P.Base_id]
  have hdeg : P.degFr ((v.inv : End A) : A ⟶ A) = 1 := by
    have h1 : P.degFr (((v.val : End A) : A ⟶ A) ≫ ((v.inv : End A) : A ⟶ A)) = 1 := by
      rw [hvv]; simp
    rw [P.degFr_comp, show P.degFr ((v.val : End A) : A ⟶ A) = 1 from hl, mul_one] at h1
    exact h1
  have hkey : P.Div ((v.inv : End A) : A ⟶ A) + P.Div ((v.val : End A) : A ⟶ A) = 0 := by
    have h2 := P.Div_comp ((v.val : End A) : A ⟶ A) ((v.inv : End A) : A ⟶ A)
    rw [hvv, P.Div_id, hbase, hdeg] at h2
    simpa [MonoidOn.map_id] using h2.symm
  exact (P.divisorial _).2 _ ⟨⟨P.Div _, P.Div _, by rw [add_comm]; exact hkey, hkey⟩, rfl⟩

/-- ★`𝒪^×(A)` の元は `𝒞^birat` でも零因子 `0`。 -/
theorem divGp_map_eq_zero_of_mem_otimes {A : C} {α : End A} (hα : α ∈ OTimes P A) :
    biratDivGp ((toBiratCat P G).map ((α : A ⟶ A))) = 0 := by
  show biratDivGp (toHomBirat (P := P) (G := G) ((α : A ⟶ A))) = 0
  rw [biratDivGp_toHomBirat, div_eq_zero_of_mem_otimes hα, toGp_zero]

/-! ## ★2. 代表元の次数と底 -/

/-- `u ∈ 𝒪^▷(A^birat)` の代表 `f` は、次数 1 で底が `c` と一致する。 -/
theorem otri_rep_base_deg (X : BiratCat P G) (u : OTri (biratPre P G) X)
    (W : IdxBirat P G (biratDown P G X)) (f : W.unop.left.obj ⟶ biratDown P G X)
    (hWf : HomBirat.mk W f = ((u : End X) : X ⟶ X)) :
    P.degFr f = 1 ∧ P.Base f = P.Base W.unop.hom.hom := by
  haveI hcb : IsIso (P.Base W.unop.hom.hom) := W.unop.hom.property.2.2
  refine ⟨?_, ?_⟩
  · have h1 : biratDeg (HomBirat.mk W f) = 1 := by rw [hWf]; exact u.2.2
    rwa [biratDeg_mk] at h1
  · have h1 : biratBase (HomBirat.mk W f) = 𝟙 _ := by
      rw [hWf]
      show (biratPre P G).Base ((u : End X) : X ⟶ X) = 𝟙 _
      rw [show (biratPre P G).Base ((u : End X) : X ⟶ X) = (biratPre P G).Base (𝟙 X) from u.2.1,
        (biratPre P G).Base_id]
    rw [biratBase_mk, sliceBaseOf_eq] at h1
    have h2 := congrArg (fun t => P.Base W.unop.hom.hom ≫ t) h1
    simp only [← Category.assoc, IsIso.hom_inv_id, Category.id_comp] at h2
    exact h2.trans (Category.comp_id _)

/-! ## ★3. 核の決定 -/

/-- ★★★**核 ⊆ 像** —— 零因子が `0` なら `𝒪^×(A)_𝒞` から来る。

★**integral**(`toGp` 単射)で `Div f = Div c` を出し、
`Definition 1.3, (vi)` `faithfulUpToUnits` を当てる。 -/
theorem otri_eq_map_of_divGp_eq_zero (X : BiratCat P G) (u : OTri (biratPre P G) X)
    (h0 : biratDivGp ((u : End X) : X ⟶ X) = 0) :
    ∃ α ∈ OTimes P (biratDown P G X),
      ((u : End X) : X ⟶ X) = (toBiratCat P G).map ((α : _ ⟶ _)) := by
  obtain ⟨W, f, hWf⟩ := HomBirat.exists_rep ((u : End X) : X ⟶ X)
  have hc : IsCoAngular P W.unop.hom.hom := W.unop.hom.property.1
  have hcs : IsPreStep P W.unop.hom.hom := W.unop.hom.property.2
  haveI hcb : IsIso (P.Base W.unop.hom.hom) := hcs.2
  obtain ⟨hfd, hfb⟩ := otri_rep_base_deg X u W f hWf
  -- ★因子の一致(ここで integral が効く)
  have hdiv : P.Div f = P.Div W.unop.hom.hom := by
    have h1 : sliceDivGpOf (P := P) W.unop.hom.hom hcb f = 0 := by
      rw [← biratDivGp_mk W f, hWf]; exact h0
    rw [sliceDivGpOf_eq] at h1
    have h2 : toGp _ (P.Div f) - ((P.degFr f : ℕ+) : ℕ) • toGp _ (P.Div W.unop.hom.hom) = 0 := by
      refine gpMap_map_injective_of_isIso (Φ := Φ) (inv (P.Base W.unop.hom.hom)) ?_
      rw [h1, map_zero]
    rw [hfd] at h2
    simp only [PNat.one_coe, one_smul, sub_eq_zero] at h2
    exact (P.divisorial _).1.1 h2
  -- ★`f` も co-angular pre-step((iii)(b) と底の一致から)
  have hcaf : IsCoAngular P f := G.core.coAngularOfPreStep W.unop.hom.hom hc hcs f
  have hpsf : IsPreStep P f := ⟨hfd, by show IsIso (P.Base f); rw [hfb]; exact hcb⟩
  -- ★(vi) を当てる
  obtain ⟨α, hα, hfeq⟩ := G.core.faithfulUpToUnits f W.unop.hom.hom hfb hdiv hcaf hpsf hc hcs
  refine ⟨α, hα, ?_⟩
  haveI : IsIso ((toBiratCat P G).map W.unop.hom.hom) := birat_isIso_of_coaPre _ hc hcs
  refine (cancel_epi ((toBiratCat P G).map W.unop.hom.hom)).mp ?_
  have hcu : (toBiratCat P G).map W.unop.hom.hom ≫ ((u : End X) : X ⟶ X)
      = (toBiratCat P G).map f := by
    rw [← hWf]; exact birat_toHom_comp_mk W.unop.hom.hom hc hcs f
  rw [hcu, hfeq, Functor.map_comp]
  rfl

/-- ★★★**核の決定**(葉 `ker-eq-image`)——
`biratDivGp u = 0` ⟺ `u` は `𝒪^×(A)_𝒞` の像。

★★これで `𝒪^×(A^birat)` は
`1 → 𝒪^×(A)_𝒞 の像 → 𝒪^×(A^birat) → Φ^gp(A)` という**中心拡大**の形になった
(核が中心に入ることは 2026-08-17 に取得済み)。 -/
theorem otriDivGp_eq_zero_iff (X : BiratCat P G) (u : OTri (biratPre P G) X) :
    biratDivGp ((u : End X) : X ⟶ X) = 0 ↔
      ∃ α ∈ OTimes P (biratDown P G X),
        ((u : End X) : X ⟶ X) = (toBiratCat P G).map ((α : _ ⟶ _)) := by
  refine ⟨otri_eq_map_of_divGp_eq_zero X u, ?_⟩
  rintro ⟨α, hα, hu⟩
  rw [hu]
  exact divGp_map_eq_zero_of_mem_otimes hα

/-! ## ★★`Proposition 4.4, (iii)` の「全射 ＋ 核」

原文 (FrdI p.83):
> induces, for each Abirat ∈Ob(Cbirat), a surjection O×(Abirat)  Φbirat(Abirat),

★`phiBiratAt A` を**`𝒪^×(A^birat)` の像**として定義してあるので、
**全射性は定義そのもの**である。★核の方が内容を持つ。 -/

/-- ★★★★**[FrdI] Proposition 4.4, (iii)** の全射——**定義そのもの**。 -/
theorem mem_phiBiratAt_iff (A : BiratCat P G)
    (x : Gp (Φ.val (P.toElem.obj (biratDown P G A)).base)) :
    x ∈ phiBiratAt P G A ↔
      ∃ δ ∈ OTimes (biratPre P G) A, biratDivGp ((δ : End A) : A ⟶ A) = x := Iff.rfl

/-- ★★★★**[FrdI] Proposition 4.4, (iii)** の**核**。

原文の「whose kernel is the image, via the injection `𝒪^▷(A)^gp → 𝒪^×(A^birat)`
of (ii), of `𝒪^×(A) ⊆ 𝒪^▷(A)^gp`」そのもの。

★(ii) の単射は `otri_toBirat_injective` / `otri_gp_injective`で取ってあるので、
ここでは「像に入る」ことを直接述べる。 -/
theorem phiBiratAt_ker (A : BiratCat P G) {u : End A}
    (hu : u ∈ OTimes (biratPre P G) A) :
    biratDivGp ((u : A ⟶ A)) = 0 ↔
      ∃ α ∈ OTimes P (biratDown P G A),
        (u : A ⟶ A) = (toBiratCat P G).map ((α : _ ⟶ _)) :=
  otriDivGp_eq_zero_iff A ⟨u, OTimes_le_OTri _ _ hu⟩

/-- ★locator —— `Proposition 4.4, (iii)` の全射＋核の条。 -/
def phiBiratAt_ker.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — 𝒪^×(A^birat) ↠ Φ^birat の核",
    sectionId := "frdi-prop-4-4" }

/-! ## ★★`Φ^birat` は `𝒞^birat` の同型で移り合う

★★★これが `(iii)` の「**`𝒟` 上の**部分関手」の要である。
`phiBiratAt` は `𝒞` の対象で添字付けされているが、
原文の `Φ^birat` は `𝒟` 上の関手なので、**持ち上げ方に依らない**ことが要る。 -/

/-- ★★★★**`Φ^birat` は同型で輸送される**。

★手: `δ ∈ 𝒪^×(B)` を `e ≫ δ ≫ e⁻¹ ∈ 𝒪^×(A)` へ共役し、
`biratDivGp_comp'` を 2 回使って
`Div^gp(e ≫ δ ≫ e⁻¹) = Φ^gp(Base e)(Div^gp δ)` を出す。
★鍵は `Div^gp(𝟙) = 0` から得る `Φ^gp(Base e)(Div^gp e⁻¹) + Div^gp e = 0`。 -/
theorem phiBiratAt_map_le {A B : BiratCat P G} (e : A ≅ B) :
    (phiBiratAt P G B).map (gpMap _ (Φ.map (biratBase e.hom))) ≤ phiBiratAt P G A := by
  rintro _ ⟨_, ⟨δ, hδ, rfl⟩, rfl⟩
  have hd1 : (biratPre P G).degFr e.hom = 1 := degFr_of_isIso (biratPre P G) e.hom
  have hd2 : (biratPre P G).degFr e.inv = 1 := degFr_of_isIso (biratPre P G) e.inv
  have hbe : (biratPre P G).Base e.hom ≫ (biratPre P G).Base e.inv
      = (biratPre P G).Base (𝟙 A) := by
    rw [← (biratPre P G).Base_comp, e.hom_inv_id]
  -- ★共役
  refine ⟨e.hom ≫ ((δ : End B) : B ⟶ B) ≫ e.inv, ⟨⟨?_, ?_⟩, ?_⟩, ?_⟩
  · show (biratPre P G).Base _ = (biratPre P G).Base (𝟙 A)
    rw [(biratPre P G).Base_comp, (biratPre P G).Base_comp,
      show (biratPre P G).Base ((δ : End B) : B ⟶ B) = (biratPre P G).Base (𝟙 B) from hδ.1.1,
      (biratPre P G).Base_id, Category.id_comp]
    exact hbe
  · show (biratPre P G).degFr _ = 1
    rw [(biratPre P G).degFr_comp, (biratPre P G).degFr_comp, hd1, hd2,
      show (biratPre P G).degFr ((δ : End B) : B ⟶ B) = 1 from hδ.1.2, one_mul, mul_one]
  · haveI : IsIso ((δ : End B) : B ⟶ B) :=
      (CategoryTheory.isUnit_iff_isIso ((δ : End B) : B ⟶ B)).mp hδ.2
    exact (CategoryTheory.isUnit_iff_isIso
      (e.hom ≫ ((δ : End B) : B ⟶ B) ≫ e.inv)).mpr inferInstance
  · -- ★Div^gp の計算
    have hone : gpMap _ (Φ.map (biratBase e.hom)) (biratDivGp e.inv) + biratDivGp e.hom = 0 := by
      have h := biratDivGp_comp' e.hom e.inv
      rw [e.hom_inv_id, biratDivGp_id,
        show ((biratDeg e.inv : ℕ+) : ℕ) = 1 from by rw [show biratDeg e.inv = 1 from hd2]; rfl,
        one_smul] at h
      exact h.symm
    have h2 : biratDivGp (((δ : End B) : B ⟶ B) ≫ e.inv)
        = biratDivGp e.inv + biratDivGp ((δ : End B) : B ⟶ B) := by
      have h := biratDivGp_comp' ((δ : End B) : B ⟶ B) e.inv
      rw [gpMap_biratBase_of_baseIdentity hδ.1.1,
        show ((biratDeg e.inv : ℕ+) : ℕ) = 1 from by rw [show biratDeg e.inv = 1 from hd2]; rfl,
        one_smul] at h
      exact h
    have h3 := biratDivGp_comp' e.hom (((δ : End B) : B ⟶ B) ≫ e.inv)
    rw [h2, map_add,
      show ((biratDeg (((δ : End B) : B ⟶ B) ≫ e.inv) : ℕ+) : ℕ) = 1 from by
        rw [show biratDeg (((δ : End B) : B ⟶ B) ≫ e.inv) = 1 from by
          show (biratPre P G).degFr (((δ : End B) : B ⟶ B) ≫ e.inv) = 1
          rw [(biratPre P G).degFr_comp,
            show (biratPre P G).degFr ((δ : End B) : B ⟶ B) = 1 from hδ.1.2, hd2, one_mul]]; rfl,
      one_smul] at h3
    rw [h3, add_assoc, add_left_comm, hone, add_zero]

/-- ★locator —— `Proposition 4.4, (iii)` の部分関手性(同型による輸送)。 -/
def phiBiratAt_map_le.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — Φ^birat の同型による輸送",
    sectionId := "frdi-prop-4-4" }

/-! ## ★出典の紐付け(条つき) -/

def otriDivGp_eq_zero_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 85, item := "Proposition 4.4, (ii)",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
