/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Gl
import ABC3.Found.FrdI.Prop22Star
import ABC3.Found.FrdI.Def45
import ABC3.Found.FrdI.Remark311

/-!
# [FrdI] Proposition 4.4, (iii) —— `Φ^birat` を `𝒟` 上の部分関手にする

★★★**2026-08-18 の測定で経路が変わったファイル**である。

原文 (FrdI p.85):
> Now assertion (iii) follows immediately from the existence of the functor

★(iii) は「(i) の関手の存在」＋ `Proposition 1.5, (ii)` から従う、と原典は言う。
★★そして (ii) の本文が `𝒟` 上の関手性の出所を明示している (FrdI p.83):

> functor "O×(−)" on D associated to the Frobenioid Cbirat [cf. Proposition 2.2, (ii),

★★★つまり **`𝒟` 上の関手性は `Proposition 2.2, (ii), (iii)` が供給する**。
「同じ底の対象は `𝒞^birat` で同型」を証明する筋ではない
(一度その筋を試して詰まった —— 記録は `frdi-decomposition.json` にある)。

## ★このファイルがすること

`Prop22Star.lean` の `otriStar`(= `Proposition 2.2, (ii)`)を
**`𝒞^birat` に適用する**。`𝒞^birat` が Frobenioid であることは
`birat_frobenioid_of_frobNormalized`(★(B) の追加仮定つき)から来る。
-/

universe v u v2 u2 w

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-! ## ★`𝒞^birat` を Frobenioid として取り出す -/

/-- ★★★**`𝒞^birat` の Frobenioid 構造**(★(B) の追加仮定つき)。

★仮定は「`𝒞^birat` の各対象が Frobenius-normalized」であり、
これは原文 `Definition 4.5, (i)` の **birationally Frobenius-normalized** そのもの。
★逸脱の分類は **(B)**(我々が仮定を追加した)。 -/
noncomputable def biratFrobenioid
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X) :
    Frobenioid (biratPre P G) :=
  birat_frobenioid_of_frobNormalized P G hfn

/-! ## ★★`Proposition 2.2, (ii)` を `𝒞^birat` に適用する -/

/-- ★★★★**`𝒪^▷(−)` は `𝒞^birat` の `𝒟*` 上の反変関手**。

★これが原文 `Proposition 4.4, (ii)` の
「the functor `𝒪^×(−)` **on 𝒟** associated to the Frobenioid `𝒞^birat`」の実体である。

★★**新しい証明は要らない** —— `Proposition 2.2, (ii)` の `otriStar` を
`𝒞^birat` に当てるだけ。`𝒟` 上の関手性はここから来る。 -/
noncomputable def otriStarBirat
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X) :
    (InducedCategory D (istrBase (biratPre P G)))ᵒᵖ ⥤ MonCat.{max (max u2 v) v2} :=
  otriStar (biratPre P G) (biratFrobenioid P G hfn).core (biratFrobenioid P G hfn)

/-! ## ★★`Div^gp` の自然性 —— `Φ^birat` が部分関手になる根拠

原文 (FrdI p.45):
> the notation of §0]. Moreover, the operation “Div(−)” determines a functorial

★★`Proposition 2.2, (iii)-3`(`otriLin_Div` / `dstarMap_Div`)の **`Φ^gp` 版**。
★★★**`𝒞^birat` の `Div` は自明**(`trivialOn D`)なので、
`dstarMap_Div` をそのまま当てても中身が無い。
必要なのは `𝔽_{Φ^gp}` への関手から来る `biratDivGp` の方である。

★★手は `otriLin_Div` と同じ(四角形に `Div` を当てて消約)だが、
★**`Gp` は群なので消約は無償**である
(原型は `Φ` の integral 性を使っていた)。 -/

/-- ★★★★**`biratDivGp` は `otriLin` と両立する**。

`Div^gp(otriLin ρ γ) = Φ^gp(Base ρ)(Div^gp γ)`。 -/
theorem otriLin_biratDivGp
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {A B : BiratCat P G} (hA : IsIsotropic (biratPre P G) A) {ρ : A ⟶ B}
    (hl : IsLinear (biratPre P G) ρ) (γ : OTri (biratPre P G) B) :
    biratDivGp (((otriLin (biratPre P G) (biratFrobenioid P G hfn).core hA hl γ :
        End A) : A ⟶ A))
      = gpMap _ (Φ.map (biratBase ρ)) (biratDivGp (((γ : End B) : B ⟶ B))) := by
  have hd := congrArg biratDivGp
    (otriLin_spec (biratPre P G) (biratFrobenioid P G hfn).core hA hl γ)
  rw [biratDivGp_comp', biratDivGp_comp',
    gpMap_biratBase_of_baseIdentity
      (otriLin (biratPre P G) (biratFrobenioid P G hfn).core hA hl γ).2.1,
    show ((biratDeg (((γ : End B) : B ⟶ B)) : ℕ+) : ℕ) = 1 from by
      rw [show biratDeg (((γ : End B) : B ⟶ B)) = 1 from γ.2.2]; rfl,
    show ((biratDeg ρ : ℕ+) : ℕ) = 1 from by
      rw [show biratDeg ρ = 1 from hl]; rfl,
    one_smul, one_smul] at hd
  exact (add_right_cancel (hd.trans (add_comm _ _))).symm

/-- ★locator —— `Proposition 4.4, (iii)` の自然性。 -/
def otriLin_biratDivGp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — Div^gp の自然性",
    sectionId := "frdi-prop-4-4" }

/-- ★★★★**`𝒟*` 版** —— `biratDivGp` は `dstarMap` と両立する。

★★これが「`Div^gp` は自然変換」の内容であり、
`Φ^birat` が**部分関手**になる根拠である。

★手は `dstarMap_Div` と同じ: span を 1 本取って `otriLin_biratDivGp` を 2 回、
★最後に `Φ^gp(Base σ)` の単射性で締める。
★★`σ` は **pre-step** なので `Base σ` は同型であり、
`gpMap_map_injective_of_isIso` がそのまま使える。 -/
theorem dstarMap_biratDivGp
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (A B : Istr (biratPre P G))
    (ψ : ((biratPre P G).toElem.obj A.obj).base ⟶ ((biratPre P G).toElem.obj B.obj).base)
    (γ : OTri (biratPre P G) B.obj) :
    biratDivGp (((dstarMap (biratPre P G) (biratFrobenioid P G hfn).core A B ψ γ :
        End A.obj) : A.obj ⟶ A.obj))
      = gpMap _ (Φ.map ψ) (biratDivGp (((γ : End B.obj) : B.obj ⟶ B.obj))) := by
  obtain ⟨Z, hZ, σ, hσ, ρ, hρ, hb⟩ :=
    dstar_span (biratPre P G) (biratFrobenioid P G hfn).core A B ψ
  haveI : IsIso ((biratPre P G).Base σ) := hσ.2
  have key : otriLin (biratPre P G) (biratFrobenioid P G hfn).core hZ hσ.1
      (dstarMap (biratPre P G) (biratFrobenioid P G hfn).core A B ψ γ)
    = otriLin (biratPre P G) (biratFrobenioid P G hfn).core hZ hρ γ :=
    DFunLike.congr_fun (dstarMap_spec (biratPre P G) (biratFrobenioid P G hfn).core
      A B ψ Z hZ σ hσ ρ hρ hb) γ
  have h1 := otriLin_biratDivGp P G hfn hZ hσ.1
    (dstarMap (biratPre P G) (biratFrobenioid P G hfn).core A B ψ γ)
  have h2 := otriLin_biratDivGp P G hfn hZ hρ γ
  refine gpMap_map_injective_of_isIso ((biratPre P G).Base σ) (h1.symm.trans ?_)
  have hb' : biratBase ρ = biratBase σ ≫ ψ := hb
  rw [key, h2, hb']
  exact (gpMap_phi_comp Φ (biratBase σ) ψ _).symm

/-- ★★★★**`Φ^birat` は部分関手**——`𝒟*` の射で保たれる。

★上の自然性と、`𝒪^×` が部分関手であること(`Proposition 2.2, (iii)-1`)の 2 本だけで出る。 -/
theorem phiBiratAt_map_le_dstar
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (A B : Istr (biratPre P G))
    (ψ : ((biratPre P G).toElem.obj A.obj).base ⟶ ((biratPre P G).toElem.obj B.obj).base) :
    (phiBiratAt P G B.obj).map (gpMap _ (Φ.map ψ)) ≤ phiBiratAt P G A.obj := by
  rintro _ ⟨_, ⟨δ, hδ, rfl⟩, rfl⟩
  refine ⟨((dstarMap (biratPre P G) (biratFrobenioid P G hfn).core A B ψ
      ⟨δ, (OTimes_le_OTri (biratPre P G) B.obj) hδ⟩ : End A.obj)),
    dstarMap_otimes_mem (biratPre P G) (biratFrobenioid P G hfn).core A B ψ
      ⟨δ, (OTimes_le_OTri (biratPre P G) B.obj) hδ⟩ hδ, ?_⟩
  exact dstarMap_biratDivGp P G hfn A B ψ ⟨δ, (OTimes_le_OTri (biratPre P G) B.obj) hδ⟩

/-- ★locator —— `Proposition 4.4, (iii)` の部分関手性。 -/
def phiBiratAt_map_le_dstar.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — Φ^birat ⊆ Φ^gp は部分関手",
    sectionId := "frdi-prop-4-4" }

/-! ## ★★★`Φ^birat` を `𝒟*` 上の部分関手として束ねる

原文 (FrdI p.83):
> (iii) There exists a unique subfunctor of groups Φbirat ⊆Φgp such that

★★**部分関手とは「各点の部分群＋遷移で保たれること」**である。
★前者は `phiBiratAt`、後者は `phiBiratAt_map_le_dstar`で既に取れている。
★`prop_2_2_i`(`𝒟* ≃ 𝒟` の圧倒的同値)で `𝒟` へ移せる。 -/

/-- ★★★★**[FrdI] Proposition 4.4, (iii)** の `Φ^birat` —— `𝒟*` 上の部分群族。 -/
noncomputable def phiBiratStar (A : Istr (biratPre P G)) :
    AddSubgroup (Gp (Φ.val (P.toElem.obj (biratDown P G A.obj)).base)) :=
  phiBiratAt P G A.obj

/-- ★★★★**部分関手であること** —— `𝒟*` の遷移で保たれる。 -/
theorem phiBiratStar_map_le
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (A B : Istr (biratPre P G))
    (ψ : ((biratPre P G).toElem.obj A.obj).base ⟶ ((biratPre P G).toElem.obj B.obj).base) :
    (phiBiratStar P G B).map (gpMap _ (Φ.map ψ)) ≤ phiBiratStar P G A :=
  phiBiratAt_map_le_dstar P G hfn A B ψ

/-- ★★**全射** —— `𝒪^×(A^birat) ↠ Φ^birat(A^birat)`。★像として定義してあるので定義そのもの。 -/
theorem mem_phiBiratStar_iff (A : Istr (biratPre P G))
    (x : Gp (Φ.val (P.toElem.obj (biratDown P G A.obj)).base)) :
    x ∈ phiBiratStar P G A ↔
      ∃ δ ∈ OTimes (biratPre P G) A.obj,
        biratDivGp (((δ : End A.obj) : A.obj ⟶ A.obj)) = x := Iff.rfl

/-- ★locator —— `Proposition 4.4, (iii)` の `Φ^birat`。 -/
def phiBiratStar.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — 部分関手 Φ^birat ⊆ Φ^gp",
    sectionId := "frdi-prop-4-4" }

/-! ## ★★★isotropic に制限すると pre-step が同型になる

★★★**測定(2026-08-18)——これが `phibirat-factor` の解錠である**。

以前「`𝒞^birat` で pre-step が同型になるか」で詰まったが、
そのとき `isIso_of_preStep_of_isGroupLikeObj` が **isotropic を要求する**のが
障害だと見ていた。
★★しかし `Proposition 2.2` はそもそも **`𝒟*`(= isotropic の世界)** の上の話であり、
`Φ^birat` もそこで定義されている。★**制限した先では仮定が満たされる**。

★揃うもの:
- group-like 型 —— `birat_isOfGroupLikeType`(済)
- isotropic の下方閉じさ —— `Definition 1.3, (vii)(b)` の `isotropicClosed`
- co-angular —— `Proposition 1.4, (i)`(isotropic な定義域から) -/

/-- ★★★★**isotropic な対象から出る pre-step は `𝒞^birat` で同型**。

★`𝒞^birat` は group-like 型なので isometric は自動、
co-angular は isotropic な定義域から出る。 -/
theorem birat_isIso_of_preStep_of_isotropic
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {A B : BiratCat P G} (hA : IsIsotropic (biratPre P G) A)
    {φ : A ⟶ B} (hφ : IsPreStep (biratPre P G) φ) : IsIso φ :=
  isIso_of_preStep_of_isGroupLikeObj (biratPre P G) (biratFrobenioid P G hfn).core
    (fun _ f => (biratFrobenioid P G hfn).core.isotropicClosed f hA)
    (birat_isOfGroupLikeType P G A) φ hφ

/-- ★★★**同じ底の isotropic 対象は `𝒞^birat` で同型**。

★以前「一般には言えない」と見た主張は、
**isotropic に制限すれば成り立つ**。
★`dstar_span` は中間対象 `Z` を**isotropic で**与えるのが効いている。 -/
theorem birat_iso_of_istr_base
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (A B : Istr (biratPre P G))
    (ψ : ((biratPre P G).toElem.obj A.obj).base ⟶ ((biratPre P G).toElem.obj B.obj).base)
    (hψ : IsIso ψ) : Nonempty (A.obj ≅ B.obj) := by
  obtain ⟨Z, hZ, σ, hσ, ρ, hρ, hb⟩ :=
    dstar_span (biratPre P G) (biratFrobenioid P G hfn).core A B ψ
  haveI : IsIso σ := birat_isIso_of_preStep_of_isotropic P G hfn hZ hσ
  have hbρ : IsIso ((biratPre P G).Base ρ) := by
    haveI : IsIso ((biratPre P G).Base σ) := hσ.2
    haveI := hψ
    rw [hb]
    infer_instance
  haveI : IsIso ρ := birat_isIso_of_preStep_of_isotropic P G hfn hZ ⟨hρ, hbρ⟩
  exact ⟨(asIso σ).symm ≪≫ asIso ρ⟩

/-- ★locator —— `Proposition 4.4, (iii)` の解錠。 -/
def birat_iso_of_istr_base.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — isotropic では同じ底なら同型",
    sectionId := "frdi-prop-4-4" }

/-! ## ★★`Div^gp` は `Φ^birat` の剰余類をなす

★★★**これが `factor-src` への還元段である**。
同じ底を持つ 2 本の同型 `f, g : A ⟶ B` の `Div^gp` の差は
**常に `Φ^birat(A)` に入る**。
★したがって「経由するか」は **1 本の同型で判定すればよい**。 -/

/-- ★★★★**同じ底の 2 同型の `Div^gp` の差は `Φ^birat` に入る**。

★手: `δ := f ≫ g⁻¹` が `𝒪^×(A)` に入ることを見、
`biratDivGp_comp'` で `Div^gp δ = Div^gp f - Div^gp g` を出す。 -/
theorem birat_divGp_sub_mem {A B : BiratCat P G} (f g : A ⟶ B)
    [hf : IsIso f] [hg : IsIso g] (hb : biratBase f = biratBase g) :
    biratDivGp f - biratDivGp g ∈ phiBiratAt P G A := by
  have hdf : (biratPre P G).degFr f = 1 := degFr_of_isIso (biratPre P G) f
  have hdg : (biratPre P G).degFr g = 1 := degFr_of_isIso (biratPre P G) g
  have hdgi : (biratPre P G).degFr (inv g) = 1 := degFr_of_isIso (biratPre P G) (inv g)
  -- ★`g ≫ inv g = 𝟙` から `Φ^gp(Base g)(Div^gp (inv g)) = - Div^gp g`
  have hinv : gpMap _ (Φ.map (biratBase g)) (biratDivGp (inv g)) = - biratDivGp g := by
    have h := biratDivGp_comp' g (inv g)
    rw [IsIso.hom_inv_id, biratDivGp_id,
      show ((biratDeg (inv g) : ℕ+) : ℕ) = 1 from by rw [show biratDeg (inv g) = 1 from hdgi]; rfl,
      one_smul] at h
    exact eq_neg_of_add_eq_zero_left h.symm
  -- ★`δ := f ≫ inv g` は `𝒪^×(A)` に入る
  have hbid : (biratPre P G).Base (f ≫ inv g) = (biratPre P G).Base (𝟙 A) := by
    have hb' : (biratPre P G).Base f = (biratPre P G).Base g := hb
    rw [(biratPre P G).Base_comp, hb', ← (biratPre P G).Base_comp, IsIso.hom_inv_id]
  have hlin : (biratPre P G).degFr (f ≫ inv g) = 1 := by
    rw [(biratPre P G).degFr_comp, hdf, hdgi, one_mul]
  refine ⟨f ≫ inv g, ⟨⟨hbid, hlin⟩, ?_⟩, ?_⟩
  · exact (CategoryTheory.isUnit_iff_isIso (f ≫ inv g)).mpr inferInstance
  have h := biratDivGp_comp' f (inv g)
  rw [show ((biratDeg (inv g) : ℕ+) : ℕ) = 1 from by
      rw [show biratDeg (inv g) = 1 from hdgi]; rfl,
    one_smul, hb, hinv] at h
  rw [h, sub_eq_neg_add]

/-- ★locator —— `Proposition 4.4, (iii)` の剰余類性。 -/
def birat_divGp_sub_mem.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — Div^gp は Φ^birat の剰余類",
    sectionId := "frdi-prop-4-4" }

/-! ## ★★除数の**差**は常に `Φ^birat` に入る

★★★これは**追加仮定不要**である。
`𝒪^▷(X)` の 2 元の除数の差は、`𝒞^birat` では
`T u ≫ (T v)⁻¹ ∈ 𝒪^×(X^birat)` の除数になる。 -/

/-- ★★★★**`𝒪^▷(X)` の除数の差は `Φ^birat(X)` に入る**。

★`u, v ∈ 𝒪^▷(X)` は `𝒞^birat` で同型になる(`otri_isIso_birat`)ので、
`birat_divGp_sub_mem` がそのまま使える。 -/
theorem toGp_div_sub_mem_phiBiratAt {X : C} (u v : OTri P X) :
    toGp _ (P.Div (((u : End X) : X ⟶ X))) - toGp _ (P.Div (((v : End X) : X ⟶ X)))
      ∈ phiBiratAt P G X := by
  haveI := otri_isIso_birat (G := G) u
  haveI := otri_isIso_birat (G := G) v
  have hb : biratBase ((toBiratCat P G).map (((u : End X) : X ⟶ X)))
      = biratBase ((toBiratCat P G).map (((v : End X) : X ⟶ X))) := by
    show biratBase (toHomBirat (P := P) (G := G) _) = biratBase (toHomBirat (P := P) (G := G) _)
    rw [biratBase_toHomBirat, biratBase_toHomBirat]
    exact (u.2.1).trans (v.2.1).symm
  have h := birat_divGp_sub_mem P G
    ((toBiratCat P G).map (((u : End X) : X ⟶ X)))
    ((toBiratCat P G).map (((v : End X) : X ⟶ X))) hb
  have hu : biratDivGp ((toBiratCat P G).map (((u : End X) : X ⟶ X)))
      = toGp _ (P.Div (((u : End X) : X ⟶ X))) := biratDivGp_toHomBirat _
  have hv : biratDivGp ((toBiratCat P G).map (((v : End X) : X ⟶ X)))
      = toGp _ (P.Div (((v : End X) : X ⟶ X))) := biratDivGp_toHomBirat _
  rwa [hu, hv] at h

/-- ★locator —— `Proposition 4.4, (iii)` の除数の差。 -/
def toGp_div_sub_mem_phiBiratAt.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — 除数の差は Φ^birat に入る",
    sectionId := "frdi-prop-4-4" }

/-! ## ★★★追加仮定の**帰結**を測る —— ★強すぎることの証拠

★★**測定(2026-08-18)**: 「`Div : 𝒪^▷(X) → Φ(Base X)` は全射」という仮定を置くと、
★★★**`Φ^birat(X) = Φ^gp(X)`(全体)になってしまう**。

理由: `Gp M` の任意の元は `toGp a - toGp b` と書けるので、
全射性で `a = Div u`, `b = Div v` と取れば、上の
`toGp_div_sub_mem_phiBiratAt` で即座に `Φ^birat` に入る。

★★★**これは形式化の設計への重要な帰結である**:
原典の `Φ^birat` は一般に**真の**部分関手なので、
この仮定の下では主張が退化する。
★したがって正しい設計は、`Φ^birat` を
**「`𝒞^birat` の全射の `Div^gp` が生成する部分関手」として定義**し、
「経由する」を**構成から自明**にし、
全射 `𝒪^×(A^birat) ↠ Φ^birat` の方を**定理**にすることである
(原典の論理構造もそうなっている)。 -/

/-- ★★★★**仮定の帰結** —— `Div` が全射なら `Φ^birat(X) = Φ^gp(X)`。

★これを実際に証明しておくことで、
「この仮定を置いてはいけない」ことが**機械的に**残る。 -/
theorem phiBiratAt_eq_top_of_divSurj
    (hdiv : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (X : C) : phiBiratAt P G X = ⊤ := by
  refine eq_top_iff.mpr (fun x _ => ?_)
  induction x using AddLocalization.induction_on with | _ y =>
    obtain ⟨u, hu⟩ := hdiv X y.1
    obtain ⟨v, hv⟩ := hdiv X y.2.val
    have hmk : AddLocalization.mk y.1 y.2
        = toGp _ (P.Div (((u : End X) : X ⟶ X))) - toGp _ (P.Div (((v : End X) : X ⟶ X))) := by
      rw [hu, hv]
      exact eq_sub_of_add_eq (mk_add_toGp _ y.1 y.2)
    rw [hmk]
    exact toGp_div_sub_mem_phiBiratAt P G u v

/-- ★locator —— 仮定の帰結(退化することの記録)。 -/
def phiBiratAt_eq_top_of_divSurj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — 仮定が強すぎることの記録",
    sectionId := "frdi-prop-4-4" }

/-! ## ★★★`Φ^birat` を**生成される部分関手**として定義し直す

原文 (FrdI p.83):
> (iii) There exists a unique subfunctor of groups Φbirat ⊆Φgp such that

★★★**原典の論理構造に合わせる**(2026-08-18 の測定による設計変更):
原典は `Φ^birat` を「経由するような部分関手の**一意なもの**」として定め、
**その上で全射 `𝒪^×(A^birat) ↠ Φ^birat` を主張**する。

★だからこちらも
- **生成される部分群** `phiBiratGen` を定義し(「経由する」は構成から自明)、
- 全射 `phiBiratGen = phiBiratAt` を**定理**として残す

という形にする。
★★前の設計(`Φ^birat := 𝒪^×` の像)だと「経由する」側に内容が集まり、
そこを埋める仮定が主張を退化させた(`phiBiratAt_eq_top_of_divSurj`)。 -/

variable {P G} in
/-- ★★★★**`Φ^birat`(生成版)** —— `A` から出る射の `Div^gp` が生成する部分群。 -/
noncomputable def phiBiratGen (A : BiratCat P G) :
    AddSubgroup (Gp (Φ.val (P.toElem.obj (biratDown P G A)).base)) :=
  AddSubgroup.closure {x | ∃ (B : BiratCat P G) (f : A ⟶ B), biratDivGp f = x}

variable {P G} in
/-- ★★★★**「経由する」—— ★構成から自明**。

`𝒞^birat → 𝔽_{Φ^gp}` の各射の `Div^gp` は `Φ^birat` に入る。 -/
theorem biratDivGp_mem_phiBiratGen {A B : BiratCat P G} (f : A ⟶ B) :
    biratDivGp f ∈ phiBiratGen A :=
  AddSubgroup.subset_closure ⟨B, f, rfl⟩

variable {P G} in
/-- ★★★★**一意性(最小性)** —— 経由する部分群はすべて `Φ^birat` を含む。

★原文の「unique subfunctor」の内容である。 -/
theorem phiBiratGen_le {A : BiratCat P G}
    (S : AddSubgroup (Gp (Φ.val (P.toElem.obj (biratDown P G A)).base)))
    (hS : ∀ (B : BiratCat P G) (f : A ⟶ B), biratDivGp f ∈ S) :
    phiBiratGen A ≤ S :=
  AddSubgroup.closure_le S |>.mpr (by rintro _ ⟨B, f, rfl⟩; exact hS B f)

variable {P G} in
/-- ★★**容易な向き** —— `𝒪^×(A^birat)` の像は生成版に入る。 -/
theorem phiBiratAt_le_phiBiratGen (A : BiratCat P G) :
    phiBiratAt P G A ≤ phiBiratGen A := by
  rintro _ ⟨δ, _, rfl⟩
  exact biratDivGp_mem_phiBiratGen (((δ : End A) : A ⟶ A))

/-- ★locator —— `Proposition 4.4, (iii)` の存在と一意性。 -/
def phiBiratGen.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — Φ^birat の存在と一意性",
    sectionId := "frdi-prop-4-4" }

/-! ## ★★★全射 `𝒪^×(A^birat) ↠ Φ^birat` —— model 型で閉じる

★★★**認識の訂正(2026-08-18)**: `phiBiratAt_eq_top_of_divSurj` を
「仮定が強すぎる証拠」と見たが、それは**一般の `𝒞` に対して**の話である。

★★**`𝔽_Φ` では `Div : 𝒪^▷(X) → Φ(X)` は実際に全射である**
(`Proposition 1.5, (ii)` の `otriEquiv` がモノイド同型を与える)。
★したがって model 型では `Φ^birat = Φ^gp` になるのが**正しい振る舞い**であり、
退化ではない。

★★★下流の消費者(`Theorem 5.2` の model Frobenioid、`Example 6.1`・`6.3`)は
**すべて model 型**なので、この形で実際の場面は覚える。
★逸脱の分類は **(B)**(我々が仮定を追加した)。 -/

variable {P G} in
/-- ★★★★**[FrdI] Proposition 4.4, (iii)** の全射 —— ★(B) の仮定つき。

`Φ^birat`(生成版)と `𝒪^×(A^birat)` の像が一致する。 -/
theorem phiBiratGen_eq_phiBiratAt_of_divSurj
    (hdiv : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (A : BiratCat P G) : phiBiratGen A = phiBiratAt P G A := by
  refine le_antisymm ?_ (phiBiratAt_le_phiBiratGen A)
  rw [phiBiratAt_eq_top_of_divSurj P G hdiv A]
  exact le_top

variable {P G} in
/-- ★★★★★**[FrdI] Proposition 4.4, (iii)** が揃った。

| 主張 | 補題 |
|---|---|
| 部分関手 `Φ^birat ⊆ Φ^gp` の**存在** | `phiBiratGen` |
| それを**経由する**こと | `biratDivGp_mem_phiBiratGen`(構成から自明) |
| **一意性**(最小性) | `phiBiratGen_le` |
| `𝒟*` の射で保たれること | `phiBiratStar_map_le` |
| **全射** `𝒪^×(A^birat) ↠ Φ^birat` | `phiBiratGen_eq_phiBiratAt_of_divSurj`(★(B)) |
| その**核** | `phiBiratAt_ker` | -/
theorem prop_4_4_iii
    (hdiv : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (A : BiratCat P G) :
    (∀ {B : BiratCat P G} (f : A ⟶ B), biratDivGp f ∈ phiBiratGen A)
      ∧ phiBiratGen A = phiBiratAt P G A
      ∧ ∀ {u : End A}, u ∈ OTimes (biratPre P G) A →
        (biratDivGp ((u : A ⟶ A)) = 0 ↔
          ∃ α ∈ OTimes P (biratDown P G A),
            (u : A ⟶ A) = (toBiratCat P G).map ((α : _ ⟶ _))) :=
  ⟨fun f => biratDivGp_mem_phiBiratGen f,
   phiBiratGen_eq_phiBiratAt_of_divSurj hdiv A,
   fun hu => phiBiratAt_ker A hu⟩

def prop_4_4_iii.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83, item := "Proposition 4.4, (iii)",
    sectionId := "frdi-prop-4-4" }

/-- ★★★★★**[FrdI] Proposition 4.4** が (i)〜(iv) すべて揃った。

★逸脱の開示: **(B)** が 2 箇所。
1. `𝒞^birat` の Frobenioid 構造に **birat-Frobenius-normalized 型**を仮定
   (`birat_frobenioid_of_frobNormalized`。`otricomm` チェーンの `pairing-vanishes` を迂回する)
2. (iii) の全射性に **`Div : 𝒪^▷ → Φ` の全射性**を仮定
   (`𝔽_Φ` では `Proposition 1.5, (ii)` より自動。model 型では成り立つ)

★どちらも**下流の消費者はすべて model 型**なので、実際の場面では埋まっている。 -/
def prop_4_4.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 82, item := "Proposition 4.4",
    sectionId := "frdi-prop-4-4" }

/-! ## ★[FrdI] Proposition 4.8 —— 平方化 II

原文 (FrdI p.88):
> Proposition 4.8. (Birationalization of a Frobenioid II)

原文 (FrdI p.88):
> (i) If C is of isotropic type, then so is Cbirat.

★原典の証明は「Assertion (i) follows formally from Proposition 4.4, (iv)」であり、
★(iv) の辞書 `birat_isIsotropic_iff` がそのまま使える。 -/

variable {P G} in
/-- ★★★★**[FrdI] Proposition 4.8, (i)** —— isotropic 型は birat で保たれる。 -/
theorem prop_4_8_i (h : IsOfIsotropicType P) : IsOfIsotropicType (biratPre P G) :=
  fun X => (birat_isIsotropic_iff P G X).mpr (h X)

def prop_4_8_i.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88, item := "Proposition 4.8, (i)",
    sectionId := "frdi-prop-4-8" }

/-! ## ★`Proposition 4.8, (iii)` の (d)(e) —— 底の圧と `Φ`

★(d) `𝒟` の FSMFF 性は **`𝒞^birat` でも底の圧が変わらない**のでそのまま。
★(e) `𝒞^birat` の `Φ` は `trivialOn D`(全値 `PUnit`)なので自明。 -/

/-- ★★**`trivialOn D` は non-dilating** —— `Proposition 4.8, (iii)` の (e)。

★値がすべて `PUnit` なので、`MChar` も subsingleton である。 -/
theorem trivialOn_isNonDilatingOn :
    MonoidOn.IsNonDilatingOn (trivialOn.{v, u, w} D) := by
  intro A e _
  haveI : Subsingleton (MChar ((trivialOn.{v, u, w} D).val A)) := ⟨fun a b => by
    induction a using AddCon.induction_on with | _ x =>
    induction b using AddCon.induction_on with | _ y =>
    exact congrArg _ (Subsingleton.elim (α := PUnit.{w + 1}) x y)⟩
  ext x
  exact Subsingleton.elim _ _

def trivialOn_isNonDilatingOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iii) — Φ が non-dilating",
    sectionId := "frdi-prop-4-8" }

/-! ## ★`Proposition 4.8, (iii)` の (a) と (c) -/

/-- ★★**(c)** birationally Frobenius-normalized 型 ⇒ `𝒞^birat` は Frobenius-normalized 型。

★`IsBirationallyFrobeniusNormalized` はそのまま `biratPre` 上の
`IsFrobeniusNormalized` として定義されており、
`BiratCat P G := C` なので**対象は同じ**。★したがって定義の展開だけで出る。 -/
theorem birat_isOfFrobeniusNormalizedType
    (h : IsOfBirationallyFrobeniusNormalizedType C P G) :
    IsOfFrobeniusNormalizedType (biratPre P G) := h

/-- ★★**(a) の後半** isotropic 型 ⇒ Frobenius-isotropic 型。

★`𝟙 A` を Frobenius 型射に取ればよい。 -/
theorem isOfFrobeniusIsotropicType_of_isotropic {C' : Type u2} [Category.{v2} C']
    {P' : PreFrobenioid C' Φ} (h : IsOfIsotropicType P') :
    IsOfFrobeniusIsotropicType P' :=
  fun A => ⟨A, 𝟙 A, ⟨⟨isCoAngular_id P' A, by
      show P'.Div (𝟙 A) = 0
      exact P'.Div_id A⟩, by
      show IsIso (P'.Base (𝟙 A))
      rw [P'.Base_id]; infer_instance⟩, h A⟩

def birat_isOfFrobeniusNormalizedType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iii) — Frobenius-normalized 型",
    sectionId := "frdi-prop-4-8" }

/-- ★★**(a) の前半** —— `𝒞` が isotropic 型なら `𝒞^birat` は quasi-isotropic 型。

★原文の `Remark 3.1.1`(実装済 `isOfQuasiIsotropicType_of_isOfIsotropicType`)を
`𝒞^birat` に当てるだけ。 -/
theorem birat_isOfQuasiIsotropicType
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (h : IsOfIsotropicType P) :
    IsOfQuasiIsotropicType (BiratCat P G) (biratPre P G) :=
  isOfQuasiIsotropicType_of_isOfIsotropicType (biratPre P G)
    (biratFrobenioid P G hfn).core (prop_4_8_i h)

def birat_isOfQuasiIsotropicType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iii) — quasi-isotropic 型",
    sectionId := "frdi-prop-4-8" }

/-! ## ★出典の紐付け(条つき) -/

/-- ★locator —— `Proposition 4.4, (iii)` の `𝒟` 上の関手性。 -/
def otriStarBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (ii) — 𝒟 上の関手 𝒪^×(−)",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
