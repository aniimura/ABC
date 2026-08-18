/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Gl
import ABC3.Found.FrdI.Prop22Star

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

/-! ## ★出典の紐付け(条つき) -/

/-- ★locator —— `Proposition 4.4, (iii)` の `𝒟` 上の関手性。 -/
def otriStarBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (ii) — 𝒟 上の関手 𝒪^×(−)",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
