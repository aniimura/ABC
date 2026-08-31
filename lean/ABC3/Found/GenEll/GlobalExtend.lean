/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.OpenExtend

/-!
# [GenEll] Remark 1.5.1 —— **固有性から点を大域的に延ばす**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

## ★★★★★★★★★★段 2c-2 —— 貼り合わせ

`OpenExtend.lean` で各素点 `v` の開近傍 `U_v` と射 `g_v : U_v ⟶ X′` を得た。
★本ファイルはそれを**貼り合わせて `Spec 𝓞_F ⟶ X′` を作る**。

| 部品 | どこから |
|---|---|
| 被覆 `⨆ U_v = ⊤` | `iSup_opens_eq_top`（`OpenExtend.lean`） |
| 整合 `hf` | ★**本ファイル** ——生成点で `xK` に揃うことから |
| 貼り合わせ | mathlib の `Scheme.Cover.glueMorphisms` |
| `⊤` から `Spec 𝓞_F` へ | mathlib の `Scheme.topIso` |

## ★★★★★★整合条件の機構

`U_v ⊓ U_w`（＝引き戻し）は `Spec 𝓞_F` の開部分スキームなので**被約**であり、
`X′` は `Spec ℤ` 上**分離的**である。
★そこで mathlib の `ext_of_isDominant_of_isSeparated` が使える
——**生成点で一致すれば全体で一致する**。

★★生成点が交わりの中で**稠密**であることは、
`genPt` が稠密（`isDominant_specMap_fractionRing`）で
交わり ⟶ `Spec 𝓞_F` が**開埋め込み**であることから
`IsDominant.of_comp_of_isOpenImmersion` で出る。

★★★生成点で `g_v` と `g_w` が揃うことは `generic_comp_open_extend`
（`OpenExtend.lean`）——**どの `U_v` で見ても `xK` に戻る**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits IsDedekindDomain

variable (F : Type) [Field F] [NumberField F]

/-! ## ★★★★★★生成点は被覆の射で一致する -/

/-- ★★★★**生成点は開近傍の中で見ても、被覆の射で一致する**。

★機構は `(⨆ U).ι` が**単射（mono）**であること
——どちらも `Spec 𝓞_F` へ落とすと素直な生成点になる。 -/
theorem generic_compat (U : FinitePlace F → (specRingOfIntegers F).Opens)
    (hxU : ∀ v, placePoint F v ∈ U v) (v w : FinitePlace F) :
    (Spec.map (stalkToF F v) ≫ (U v).fromSpecStalkOfMem _ (hxU v)) ≫
        (specRingOfIntegers F).homOfLE (le_iSup U v)
      = (Spec.map (stalkToF F w) ≫ (U w).fromSpecStalkOfMem _ (hxU w)) ≫
        (specRingOfIntegers F).homOfLE (le_iSup U w) := by
  refine (cancel_mono ((⨆ z, U z).ι)).1 ?_
  simp only [Category.assoc]
  rw [Scheme.homOfLE_ι _ (le_iSup U v), Scheme.homOfLE_ι _ (le_iSup U w),
    ← Category.assoc, ← Category.assoc, generic_comp_ι F v (U v) (hxU v),
    generic_comp_ι F w (U w) (hxU w)]

/-! ## ★★★★★★★★整合条件 -/

/-- ★★★★★★★★**貼り合わせの整合条件**。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

★機構は 3 つ:
(a) 交わりは `Spec 𝓞_F` の開部分スキームなので**被約**、
(b) 生成点は交わりの中で**稠密**（`IsDominant.of_comp_of_isOpenImmersion`）、
(c) 生成点では両方 `xK` に戻る（`generic_comp_open_extend`）。
★★あとは `ext_of_isDominant_of_isSeparated` 1 本。 -/
theorem hf_agreement {X' : Scheme.{0}} (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsProper f']
    (xK : Spec (CommRingCat.of F) ⟶ X')
    (xv : (v : FinitePlace F) → Spec (CommRingCat.of (Localization.AtPrime v.asIdeal)) ⟶ X')
    (hxv : ∀ v, Spec.map (CommRingCat.ofHom
        (algebraMap (Localization.AtPrime v.asIdeal) F)) ≫ xv v = xK)
    (U : FinitePlace F → (specRingOfIntegers F).Opens)
    (hxU : ∀ v, placePoint F v ∈ U v)
    (g : (v : FinitePlace F) → (U v).toScheme ⟶ X')
    (hg : ∀ v, Spec.map (Spec.stalkIso (CommRingCat.of (NumberField.RingOfIntegers F))
        (placePoint F v)).inv ≫ xv v = (U v).fromSpecStalkOfMem _ (hxU v) ≫ g v)
    (v w : FinitePlace F) :
    pullback.fst ((specRingOfIntegers F).homOfLE (le_iSup U v))
      ((specRingOfIntegers F).homOfLE (le_iSup U w)) ≫ g v
      = pullback.snd _ _ ≫ g w := by
  haveI : IsReduced (pullback ((specRingOfIntegers F).homOfLE (le_iSup U v))
      ((specRingOfIntegers F).homOfLE (le_iSup U w))) :=
    isReduced_of_isOpenImmersion (pullback.fst _ _)
  set ι := pullback.lift _ _ (generic_compat F U hxU v w) with hι
  set p : pullback ((specRingOfIntegers F).homOfLE (le_iSup U v))
      ((specRingOfIntegers F).homOfLE (le_iSup U w)) ⟶ specRingOfIntegers F :=
    pullback.fst _ _ ≫ (specRingOfIntegers F).homOfLE (le_iSup U v) ≫ (⨆ z, U z).ι with hp
  have hc : ι ≫ p
      = Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) F)) := by
    rw [hι, hp]
    simp only [← Category.assoc]
    rw [pullback.lift_fst, Category.assoc, Scheme.homOfLE_ι _ (le_iSup U v)]
    exact generic_comp_ι F v (U v) (hxU v)
  haveI hd0 : IsDominant (Spec.map (CommRingCat.ofHom
      (algebraMap (NumberField.RingOfIntegers F) F))) :=
    isDominant_specMap_fractionRing (NumberField.RingOfIntegers F) F
  haveI : IsDominant (ι ≫ p) := hc ▸ hd0
  haveI : IsOpenImmersion p := by rw [hp]; infer_instance
  haveI : IsDominant ι := IsDominant.of_comp_of_isOpenImmersion ι p
  refine ext_of_isDominant_of_isSeparated f' (specZIsTerminal.hom_ext _ _) ι ?_
  rw [hι]
  simp only [← Category.assoc]
  rw [pullback.lift_fst, pullback.lift_snd,
    generic_comp_open_extend F v (xv v) (U v) (hxU v) (g v) (hg v) xK (hxv v),
    generic_comp_open_extend F w (xv w) (U w) (hxU w) (g w) (hg w) xK (hxv w)]

/-! ## ★★★★★★★★★★到達点 —— 大域延長の存在 -/

/-- ★★★★★★★★★★**[GenEll] Definition 1.5 —— 固有性から点が延びる**。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

`X′` が `ℤ` 上**固有**なら、`F`-点は `Spec 𝓞_F`-点へ延びる。

★★★★★これが原文の「where we recall that X is proper!」の中身であり、
`Remark 1.5.1` の点の対応 `ePt` を作る段でもある。

★組み立ては 5 段:
1. 付値環での延長（`ProperExtend.lean`、mathlib の付値判定法）
2a. 各素点 `v` で `(𝓞_F)_v`-点へ（`DedekindExtend.lean`）
2b. 大域延長の一意性（`DedekindExtend.lean`）
2c-1. 開近傍へ広げる（`OpenExtend.lean`、Stacks 0BX6）
2c-2. ★**貼り合わせ**（本ファイル）

★★`extend_unique`（`DedekindExtend.lean`）と合わせると、
**延長は存在して一意**である。 -/
theorem exists_global_extend {X' : Scheme.{0}} (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsProper f']
    (xK : Spec (CommRingCat.of F) ⟶ X') :
    ∃ xR : specRingOfIntegers F ⟶ X',
      Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) F)) ≫ xR = xK := by
  classical
  choose xv hxv using fun v : FinitePlace F => (exists_unique_extend_atPrime F f' v xK).exists
  choose U hxU g hg using fun v : FinitePlace F => exists_open_extend F f' v (xv v)
  have hsup := iSup_opens_eq_top F U hxU
  haveI : IsIso ((⊤ : (specRingOfIntegers F).Opens).ι) :=
    (Scheme.topIso (specRingOfIntegers F)).isIso_hom
  haveI hiso : IsIso ((⨆ v, U v).ι) := by rw [hsup]; infer_instance
  refine ⟨inv ((⨆ v, U v).ι) ≫ (Scheme.Opens.iSupOpenCover U).glueMorphisms g
    (fun v w => hf_agreement F f' xK xv hxv U hxU g hg v w), ?_⟩
  obtain ⟨v₀⟩ := nonempty_finitePlace F
  have hglue := Scheme.Cover.ι_glueMorphisms (Scheme.Opens.iSupOpenCover U) g
    (fun v w => hf_agreement F f' xK xv hxv U hxU g hg v w) v₀
  have key : Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) F))
      = (Spec.map (stalkToF F v₀) ≫ (U v₀).fromSpecStalkOfMem _ (hxU v₀)) ≫
        ((Scheme.Opens.iSupOpenCover U).f v₀ ≫ (⨆ z, U z).ι) := by
    show _ = (Spec.map (stalkToF F v₀) ≫ (U v₀).fromSpecStalkOfMem _ (hxU v₀)) ≫
      ((specRingOfIntegers F).homOfLE (le_iSup U v₀) ≫ (⨆ z, U z).ι)
    rw [Scheme.homOfLE_ι _ (le_iSup U v₀)]
    exact (generic_comp_ι F v₀ (U v₀) (hxU v₀)).symm
  rw [key]
  simp only [Category.assoc, IsIso.hom_inv_id_assoc]
  refine Eq.trans ?_
    (generic_comp_open_extend F v₀ (xv v₀) (U v₀) (hxU v₀) (g v₀) (hg v₀) xK (hxv v₀))
  simp only [Category.assoc]
  exact congrArg
    (fun z => Spec.map (stalkToF F v₀) ≫ (U v₀).fromSpecStalkOfMem (placePoint F v₀) (hxU v₀) ≫ z)
    hglue

/-- ★★★★★★★★★★**存在と一意性を合わせた形**。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

★★★これで `X′(F) ≃ X′(𝓞_F)`（`ℤ` 上固有な `X′` について）が型で書けた。 -/
theorem existsUnique_global_extend {X' : Scheme.{0}}
    (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsProper f']
    (xK : Spec (CommRingCat.of F) ⟶ X') :
    ∃! xR : specRingOfIntegers F ⟶ X',
      Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) F)) ≫ xR = xK := by
  obtain ⟨xR, hxR⟩ := exists_global_extend F f' xK
  exact ⟨xR, hxR, fun y hy => extend_unique F f' y xR (by rw [hy, hxR])⟩

/-! ### ★出典の紐付け(`.src`) -/

def exists_global_extend.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(固有性から点を延ばす——大域)",
    sectionId := "genell-def-1-5" }

def existsUnique_global_extend.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(固有なら X′(F) ≃ X′(𝓞_F))",
    sectionId := "genell-def-1-5" }

def exists_global_extend.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_open_extend(局所延長を開近傍へ広げる)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_open_extend") 8,
    .citation "[ABC3]" "iSup_opens_eq_top(開近傍は全体を覆う)"
      (.inProject "ABC3" "ABC3.Found.GenEll.iSup_opens_eq_top") 8,
    .citation "[mathlib]" "Scheme.Cover.glueMorphisms(整合する射の貼り合わせ)"
      (.inMathlib "AlgebraicGeometry.Scheme.Cover.glueMorphisms") 8,
    .citation "[mathlib]" "ext_of_isDominant_of_isSeparated(生成点で一致すれば全体で一致)"
      (.inMathlib "AlgebraicGeometry.ext_of_isDominant_of_isSeparated") 8,
    .citation "[mathlib]" "IsDominant.of_comp_of_isOpenImmersion(開埋め込みで稠密性を降ろす)"
      (.inMathlib "AlgebraicGeometry.IsDominant.of_comp_of_isOpenImmersion") 8,
    .implicitStep
      ("★★★★★原文は「where we recall that X is proper!」の 1 語で済ませている。" ++
       "★形式化では 5 段に分かれた: 付値環での延長 → 各素点での延長 → 一意性 → " ++
       "開近傍へ広げる → 貼り合わせ") 8,
    .implicitStep
      ("★整合条件の機構は 3 つ: 交わりは開部分スキームなので被約、" ++
       "生成点は交わりの中で稠密、生成点では両方 xK に戻る") 8 ]

end ABC3.Found.GenEll
