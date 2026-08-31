/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PullbackGen
import ABC3.Found.Arakelov.TensorMetric
import ABC3.Found.Arakelov.AMetricMonoid
import ABC3.Meta.Claim

/-!
# **算術直線束の引き戻し**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

## ★★★★★★★★★★「evident」を開く

原文は引き戻しを **evident**（自明）で畳んでいる。畳まれていた中身は

| 段 | 内容 | 場所 |
|---|---|---|
| 1 | 層の引き戻し `f^*L` | `PicSchemeDelta.lean`（`pullbackPre`） |
| 2 | 局所自明性が保たれること | `PicLTPull.lean` |
| 3 | 自明化の輸送 `pullTrivOfBase` | `PullbackTriv.lean` |
| 4 | 遷移単元が輸送と可換（同じ `W`） | `PullbackUnitEnd.lean`（§9-740） |
| 5 | 遷移単元が `W` の縮小に依らない | `PullbackGen.lean`（§9-741） |
| 6 | **計量の引き戻し** | ★★**本ファイル** |

## ★★★★★★★★構成 —— `LocalMetric.tensor` と同じ形

`V` 上の自明化 `e` に対して値を決めたい。★`f^*L` は `V` 上で自明でも
`L` は `f(V)` の上で自明とは限らないので、`p` の周りで

    チャート `c = (W, V')`：`V' ≤ V`、`V' ≤ f⁻¹W`、`p ∈ V'`、`e_W : L|_W ≅ 𝟙_`

を取り、そこで `m.h W e_W (p ≫ f)` と置いて遷移単元で `V` の自明化 `e` へ運ぶ:

    `pullChartH c e = ‖⟨輸送した自明化 → e⟩ の遷移単元 (p)‖⁻¹ · m.h W e_W (p ≫ f)`

## ★★★★★★★★チャート独立性の 2 つの向き

`TensorMetric` では開集合が 1 つ（`X` の上）だったが、
ここでは **`Y` の `W` と `X` の `V'` の 2 つ**を縮める必要がある。

| 段 | 何を動かすか | 機構 |
|---|---|---|
| `pullChartH_triv_indep` | `e_W`（同じ `W`） | §9-740 `transUnit_pullTrivOfBase` ＋ `compat` |
| `pullChartH_shrinkV` | `V'`（`X` 側） | `pullTrivOfBase_shrinkV`（`rfl`）＋ `transUnit_restrict` |
| `pullChartH_shrinkW` | `W`（`Y` 側） | §9-741 `transUnit_pullTrivOfBase_shrinkW` ＋ `restrict` |

★`pullChartH_indep` は `c.W ⊓ c'.W` と `c.V' ⊓ c'.V'` へ**両方**降ろす。

## ★★★★★★★★★★到達点 —— 特徴づけ

    `(f^* m).h V' (pullTrivOfBase f L W e_W) p = m.h W e_W (p ≫ f)`

★「輸送した自明化の上では、値は底の計量の値そのもの」。
★★これが `Classical.choice` の影を消す——値はチャートに依らないので、
この等式が**すべての** `W`, `e_W` について成り立つ。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll

/-! ## ★★★★★引き戻しの計量のチャート -/

/-- ★★★**引き戻しの計量のチャート**。

★`Y` の開集合 `W`（そこで `L` が自明）と `X` の開集合 `V'`（`V` の中、`f⁻¹W` の中、
`p` を含む）の**組**である。★★`TensorMetric` のチャートと違い、
開集合が 2 つ（底と上）あるのが引き戻しの特徴である。 -/
structure PullChart {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (V : X.Opens) (p : Spec (CommRingCat.of ℂ) ⟶ X) where
  /-- `Y` の側の開集合——ここで `L` が自明になる。 -/
  W : Y.Opens
  /-- `X` の側の開集合。 -/
  V' : X.Opens
  /-- `V` の中にある。 -/
  hV'V : V' ≤ V
  /-- `f⁻¹W` の中にある。 -/
  hV'W : V' ≤ (Opens.map f.base).obj W
  /-- `p` を含む。 -/
  hpV' : p ⁻¹ᵁ V' = ⊤
  /-- `L` の `W` 上の自明化。 -/
  eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W)

/-- ★★★**`L` が局所自明なら、`p` を含む `V` の中にチャートが取れる**。

★`Spec ℂ` は 1 点なので、`p ≫ f` の像の点を含む自明化の開集合 `W` を取り、
`V' := V ⊓ f⁻¹W` と置けばよい。 -/
theorem nonempty_pullChart {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    (hL : IsLocallyTrivial Y L) (V : X.Opens) (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (hp : p ⁻¹ᵁ V = ⊤) : Nonempty (PullChart f L V p) := by
  obtain ⟨S, hS, htriv⟩ := hL ⊤
  have hxT : (p ≫ f).base default ∈ (⊤ : Y.Opens) := trivial
  obtain ⟨W, g, hg, hW⟩ := hS ((p ≫ f).base default) hxT
  have hWT : W ≤ (⊤ : Y.Opens) := le_top
  have hgt : S.arrows (homOfLE hWT) := (Subsingleton.elim g (homOfLE hWT)) ▸ hg
  have hpW : p ⁻¹ᵁ ((Opens.map f.base).obj W) = ⊤ :=
    preimage_eq_top_of_mem p ((Opens.map f.base).obj W)
      (fun z => by rw [Subsingleton.elim z default]; exact hW)
  have hpV' : p ⁻¹ᵁ (V ⊓ (Opens.map f.base).obj W) = ⊤ := by
    show p ⁻¹ᵁ V ⊓ p ⁻¹ᵁ ((Opens.map f.base).obj W) = ⊤
    rw [hp, hpW, inf_idem]
  exact ⟨⟨W, V ⊓ (Opens.map f.base).obj W, inf_le_left, inf_le_right, hpV',
    (htriv (homOfLE hWT) hgt).some⟩⟩

/-- ★チャートを大きい開集合の上のチャートと見る。 -/
def PullChart.lift {X Y : Scheme.{0}} {f : X ⟶ Y} {L : Y.PresheafOfModules}
    {V V₂ : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X} (c : PullChart f L V p)
    (hVV₂ : V ≤ V₂) : PullChart f L V₂ p :=
  ⟨c.W, c.V', c.hV'V.trans hVV₂, c.hV'W, c.hpV', c.eW⟩

/-- ★チャートの `X` 側の開集合を小さくする。 -/
noncomputable def PullChart.shrinkV {X Y : Scheme.{0}} {f : X ⟶ Y} {L : Y.PresheafOfModules}
    {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X} (c : PullChart f L V p)
    {V'' : X.Opens} (hV''V' : V'' ≤ c.V') (hpV'' : p ⁻¹ᵁ V'' = ⊤) : PullChart f L V p :=
  ⟨c.W, V'', hV''V'.trans c.hV'V, hV''V'.trans c.hV'W, hpV'', c.eW⟩

/-- ★チャートの `Y` 側の開集合を小さくする。 -/
noncomputable def PullChart.shrinkW {X Y : Scheme.{0}} {f : X ⟶ Y} {L : Y.PresheafOfModules}
    {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X} (c : PullChart f L V p)
    {W'' : Y.Opens} (hW : W'' ≤ c.W) {V'' : X.Opens} (hV''V' : V'' ≤ c.V')
    (hV''W'' : V'' ≤ (Opens.map f.base).obj W'') (hpV'' : p ⁻¹ᵁ V'' = ⊤) :
    PullChart f L V p :=
  ⟨W'', V'', hV''V'.trans c.hV'V, hV''W'', hpV'', trivialOfLe L hW c.eW⟩

/-! ## ★★★★★★候補値とその選択独立性 -/

/-- ★★**チャートが与える基準ノルムの候補値**。

    `‖⟨輸送した自明化 → e⟩ の遷移単元 (p)‖⁻¹ · m.h W e_W (p ≫ f)` -/
noncomputable def pullChartH {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    (m : LocalMetric Y L) {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X}
    (c : PullChart f L V p)
    (e : (restrictPresheafFunctor X V).obj ((pullbackPre f).obj L)
      ≅ 𝟙_ (PresheafModulesOn X V)) : ℝ :=
  ‖evalOn p c.V' c.hpV' (transUnit ((pullbackPre f).obj L) c.V'
      (pullTrivOfBase f L c.W c.eW c.hV'W)
      (trivialOfLe ((pullbackPre f).obj L) c.hV'V e))‖⁻¹ * m.h c.W c.eW (p ≫ f)

theorem pullChartH_pos {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    (m : LocalMetric Y L) {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X}
    (c : PullChart f L V p)
    (e : (restrictPresheafFunctor X V).obj ((pullbackPre f).obj L)
      ≅ 𝟙_ (PresheafModulesOn X V)) :
    0 < pullChartH f m c e :=
  mul_pos
    (inv_pos.2 (lt_of_le_of_ne (norm_nonneg _) (Ne.symm (norm_ne_zero_iff.2
      (evalOn_ne_zero_of_isUnit p c.V' c.hpV' (isUnit_transUnit _ c.V' _ _))))))
    (m.pos c.W c.eW (p ≫ f))

/-- ★★★★**候補値は底の自明化の取り方に依らない**（同じ `W` の上で）。

★機構は §9-740 の `transUnit_pullTrivOfBase`（遷移単元は引き戻しと可換）と
`m.compat`、そして `evalOn_pullback`（`ρ` を通すと底の値になる）。 -/
theorem pullChartH_triv_indep {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    (m : LocalMetric Y L) {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X}
    {W : Y.Opens} {V' : X.Opens} (hV'V : V' ≤ V)
    (hV'W : V' ≤ (Opens.map f.base).obj W) (hpV' : p ⁻¹ᵁ V' = ⊤)
    (eW eW' : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    (e : (restrictPresheafFunctor X V).obj ((pullbackPre f).obj L)
      ≅ 𝟙_ (PresheafModulesOn X V)) :
    pullChartH f m (⟨W, V', hV'V, hV'W, hpV', eW'⟩ : PullChart f L V p) e
      = pullChartH f m (⟨W, V', hV'V, hV'W, hpV', eW⟩ : PullChart f L V p) e := by
  have hpW : (p ≫ f) ⁻¹ᵁ W = ⊤ := comp_preimage_eq_top_of_le f hV'W hpV'
  have hA' : m.h W eW' (p ≫ f)
      = m.h W eW (p ≫ f) * ‖evalOn (p ≫ f) W hpW (transUnit L W eW' eW)‖ := by
    rw [← m.compat W eW eW' (p ≫ f) hpW, mul_assoc,
      norm_evalOn_transUnit_symm L W eW eW' (p ≫ f) hpW, mul_one]
  have hane : ‖evalOn (p ≫ f) W hpW (transUnit L W eW' eW)‖ ≠ 0 :=
    norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit (p ≫ f) W hpW (isUnit_transUnit _ W eW' eW))
  have htne : ‖evalOn p V' hpV' (transUnit ((pullbackPre f).obj L) V'
      (pullTrivOfBase f L W eW hV'W)
      (trivialOfLe ((pullbackPre f).obj L) hV'V e))‖ ≠ 0 :=
    norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p V' hpV' (isUnit_transUnit _ V' _ _))
  have hfac : ‖evalOn p V' hpV' (transUnit ((pullbackPre f).obj L) V'
        (pullTrivOfBase f L W eW' hV'W)
        (trivialOfLe ((pullbackPre f).obj L) hV'V e))‖
      = ‖evalOn (p ≫ f) W hpW (transUnit L W eW' eW)‖
        * ‖evalOn p V' hpV' (transUnit ((pullbackPre f).obj L) V'
            (pullTrivOfBase f L W eW hV'W)
            (trivialOfLe ((pullbackPre f).obj L) hV'V e))‖ := by
    rw [transUnit_trans ((pullbackPre f).obj L) V' (pullTrivOfBase f L W eW' hV'W)
      (pullTrivOfBase f L W eW hV'W) (trivialOfLe ((pullbackPre f).obj L) hV'V e),
      transUnit_pullTrivOfBase f L W eW' eW hV'W, evalOn_mul, norm_mul,
      evalOn_pullback f p W hV'W hpV' hpW]
  show ‖evalOn p V' hpV' (transUnit ((pullbackPre f).obj L) V'
      (pullTrivOfBase f L W eW' hV'W)
      (trivialOfLe ((pullbackPre f).obj L) hV'V e))‖⁻¹ * m.h W eW' (p ≫ f)
    = ‖evalOn p V' hpV' (transUnit ((pullbackPre f).obj L) V'
        (pullTrivOfBase f L W eW hV'W)
        (trivialOfLe ((pullbackPre f).obj L) hV'V e))‖⁻¹ * m.h W eW (p ≫ f)
  rw [hfac, hA']
  field_simp

/-- ★★★★**候補値は `X` 側の開集合を小さくしても変わらない**。

★機構は `pullTrivOfBase_shrinkV`（`rfl`）・`trivialOfLe_trans`（`rfl`）と
`transUnit_restrict`・`evalOn_restrict`。★★底の値 `m.h W e_W (p ≫ f)` は動かない。 -/
theorem pullChartH_shrinkV {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    (m : LocalMetric Y L) {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X}
    (c : PullChart f L V p) {V'' : X.Opens} (hV''V' : V'' ≤ c.V') (hpV'' : p ⁻¹ᵁ V'' = ⊤)
    (e : (restrictPresheafFunctor X V).obj ((pullbackPre f).obj L)
      ≅ 𝟙_ (PresheafModulesOn X V)) :
    pullChartH f m (c.shrinkV hV''V' hpV'') e = pullChartH f m c e := by
  show ‖evalOn p V'' hpV'' (transUnit ((pullbackPre f).obj L) V''
      (pullTrivOfBase f L c.W c.eW (hV''V'.trans c.hV'W))
      (trivialOfLe ((pullbackPre f).obj L) (hV''V'.trans c.hV'V) e))‖⁻¹
      * m.h c.W c.eW (p ≫ f) = _
  rw [← pullTrivOfBase_shrinkV f L c.W c.eW hV''V' c.hV'W (hV''V'.trans c.hV'W),
    trivialOfLe_trans ((pullbackPre f).obj L) hV''V' c.hV'V e,
    transUnit_restrict ((pullbackPre f).obj L) hV''V'
      (pullTrivOfBase f L c.W c.eW c.hV'W)
      (trivialOfLe ((pullbackPre f).obj L) c.hV'V e),
    evalOn_restrict p hV''V' hpV'']
  rfl

/-- ★★★★★★**候補値は `Y` 側の開集合を小さくしても変わらない**。

★★これが引き戻し特有の段である——機構は §9-741 の
`transUnit_pullTrivOfBase_shrinkW`（生成切断が一致するので遷移単元が一致する）と
`m.restrict`。 -/
theorem pullChartH_shrinkW {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    (m : LocalMetric Y L) {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X}
    (c : PullChart f L V p) {W'' : Y.Opens} (hW : W'' ≤ c.W) {V'' : X.Opens}
    (hV''V' : V'' ≤ c.V') (hV''W'' : V'' ≤ (Opens.map f.base).obj W'')
    (hpV'' : p ⁻¹ᵁ V'' = ⊤)
    (e : (restrictPresheafFunctor X V).obj ((pullbackPre f).obj L)
      ≅ 𝟙_ (PresheafModulesOn X V)) :
    pullChartH f m (c.shrinkW hW hV''V' hV''W'' hpV'') e
      = pullChartH f m (c.shrinkV hV''V' hpV'') e := by
  have hpW'' : (p ≫ f) ⁻¹ᵁ W'' = ⊤ := comp_preimage_eq_top_of_le f hV''W'' hpV''
  show ‖evalOn p V'' hpV'' (transUnit ((pullbackPre f).obj L) V''
      (pullTrivOfBase f L W'' (trivialOfLe L hW c.eW) hV''W'')
      (trivialOfLe ((pullbackPre f).obj L) (hV''V'.trans c.hV'V) e))‖⁻¹
      * m.h W'' (trivialOfLe L hW c.eW) (p ≫ f) = _
  rw [transUnit_pullTrivOfBase_shrinkW f L hW c.eW hV''W'' (hV''V'.trans c.hV'W)
    (trivialOfLe ((pullbackPre f).obj L) (hV''V'.trans c.hV'V) e),
    m.restrict hW c.eW (p ≫ f) hpW'']
  rfl

/-- ★★★★★★★★**候補値はチャートの取り方に全く依らない**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★機構は「`c.W ⊓ c'.W` と `c.V' ⊓ c'.V'` へ**両方**降ろして `pullChartH_triv_indep`」。 -/
theorem pullChartH_indep {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    (m : LocalMetric Y L) {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X}
    (c c' : PullChart f L V p)
    (e : (restrictPresheafFunctor X V).obj ((pullbackPre f).obj L)
      ≅ 𝟙_ (PresheafModulesOn X V)) :
    pullChartH f m c e = pullChartH f m c' e := by
  have hpV'' : p ⁻¹ᵁ (c.V' ⊓ c'.V') = ⊤ := by
    show p ⁻¹ᵁ c.V' ⊓ p ⁻¹ᵁ c'.V' = ⊤
    rw [c.hpV', c'.hpV', inf_idem]
  have hV''W : c.V' ⊓ c'.V' ≤ (Opens.map f.base).obj (c.W ⊓ c'.W) :=
    fun z hz => ⟨c.hV'W hz.1, c'.hV'W hz.2⟩
  rw [← pullChartH_shrinkV f m c (inf_le_left : c.V' ⊓ c'.V' ≤ c.V') hpV'' e,
    ← pullChartH_shrinkV f m c' (inf_le_right : c.V' ⊓ c'.V' ≤ c'.V') hpV'' e,
    ← pullChartH_shrinkW f m c (inf_le_left : c.W ⊓ c'.W ≤ c.W)
      (inf_le_left : c.V' ⊓ c'.V' ≤ c.V') hV''W hpV'' e,
    ← pullChartH_shrinkW f m c' (inf_le_right : c.W ⊓ c'.W ≤ c'.W)
      (inf_le_right : c.V' ⊓ c'.V' ≤ c'.V') hV''W hpV'' e]
  exact pullChartH_triv_indep f m _ hV''W hpV''
    (trivialOfLe L (inf_le_right : c.W ⊓ c'.W ≤ c'.W) c'.eW)
    (trivialOfLe L (inf_le_left : c.W ⊓ c'.W ≤ c.W) c.eW) e

/-! ## ★★★★★★★★★★引き戻しの計量 -/

open scoped Classical in
/-- ★★★★★★★★★★**計量の引き戻し**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★`p` の周りで底が自明になるチャートを取り、そこで `m.h W e_W (p ≫ f)` と置いて
遷移単元で `V` の自明化 `e` へ運ぶ。★★値がチャートに依らないことは `pullChartH_indep`。 -/
noncomputable def LocalMetric.pullback {X Y : Scheme.{0}} (f : X ⟶ Y)
    {L : Y.PresheafOfModules} (hL : IsLocallyTrivial Y L) (m : LocalMetric Y L) :
    LocalMetric X ((pullbackPre f).obj L) where
  h V e p := if hc : Nonempty (PullChart f L V p) then pullChartH f m hc.some e else 1
  pos V e p := by
    by_cases hc : Nonempty (PullChart f L V p)
    · show 0 < (if hc : Nonempty (PullChart f L V p) then pullChartH f m hc.some e else 1)
      rw [dif_pos hc]
      exact pullChartH_pos f m hc.some e
    · show 0 < (if hc : Nonempty (PullChart f L V p) then pullChartH f m hc.some e else 1)
      rw [dif_neg hc]
      exact one_pos
  compat V e e' p hp := by
    have hc : Nonempty (PullChart f L V p) := nonempty_pullChart f hL V p hp
    show (if hc : Nonempty (PullChart f L V p) then pullChartH f m hc.some e' else 1) * _
      = (if hc : Nonempty (PullChart f L V p) then pullChartH f m hc.some e else 1)
    rw [dif_pos hc, dif_pos hc]
    set c := hc.some with hcdef
    have hstep : ‖evalOn p c.V' c.hpV' (transUnit ((pullbackPre f).obj L) c.V'
          (pullTrivOfBase f L c.W c.eW c.hV'W)
          (trivialOfLe ((pullbackPre f).obj L) c.hV'V e'))‖
        = ‖evalOn p c.V' c.hpV' (transUnit ((pullbackPre f).obj L) c.V'
            (pullTrivOfBase f L c.W c.eW c.hV'W)
            (trivialOfLe ((pullbackPre f).obj L) c.hV'V e))‖
          * ‖evalOn p V hp (transUnit ((pullbackPre f).obj L) V e e')‖ := by
      rw [transUnit_trans ((pullbackPre f).obj L) c.V' (pullTrivOfBase f L c.W c.eW c.hV'W)
        (trivialOfLe ((pullbackPre f).obj L) c.hV'V e)
        (trivialOfLe ((pullbackPre f).obj L) c.hV'V e'),
        transUnit_restrict ((pullbackPre f).obj L) c.hV'V e e', evalOn_mul, norm_mul,
        evalOn_restrict p c.hV'V c.hpV']
    have hne : ‖evalOn p c.V' c.hpV' (transUnit ((pullbackPre f).obj L) c.V'
        (pullTrivOfBase f L c.W c.eW c.hV'W)
        (trivialOfLe ((pullbackPre f).obj L) c.hV'V e))‖ ≠ 0 :=
      norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p c.V' c.hpV' (isUnit_transUnit _ c.V' _ _))
    have hune : ‖evalOn p V hp (transUnit ((pullbackPre f).obj L) V e e')‖ ≠ 0 :=
      norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p V hp (isUnit_transUnit _ V e e'))
    show ‖evalOn p c.V' c.hpV' (transUnit ((pullbackPre f).obj L) c.V'
          (pullTrivOfBase f L c.W c.eW c.hV'W)
          (trivialOfLe ((pullbackPre f).obj L) c.hV'V e'))‖⁻¹ * m.h c.W c.eW (p ≫ f)
        * ‖evalOn p V hp (transUnit ((pullbackPre f).obj L) V e e')‖
      = ‖evalOn p c.V' c.hpV' (transUnit ((pullbackPre f).obj L) c.V'
          (pullTrivOfBase f L c.W c.eW c.hV'W)
          (trivialOfLe ((pullbackPre f).obj L) c.hV'V e))‖⁻¹ * m.h c.W c.eW (p ≫ f)
    rw [hstep, mul_inv]
    field_simp
  restrict {V₁ V₂} h21 e p hp2 := by
    have hc2 : Nonempty (PullChart f L V₂ p) := nonempty_pullChart f hL V₂ p hp2
    have hc1 : Nonempty (PullChart f L V₁ p) :=
      nonempty_pullChart f hL V₁ p (preimage_eq_top_of_le h21 hp2)
    show (if hc : Nonempty (PullChart f L V₂ p) then
        pullChartH f m hc.some (trivialOfLe ((pullbackPre f).obj L) h21 e) else 1)
      = (if hc : Nonempty (PullChart f L V₁ p) then pullChartH f m hc.some e else 1)
    rw [dif_pos hc2, dif_pos hc1]
    have hlift : pullChartH f m (hc2.some.lift h21) e
        = pullChartH f m hc2.some (trivialOfLe ((pullbackPre f).obj L) h21 e) := rfl
    rw [← hlift]
    exact pullChartH_indep f m _ _ e

open scoped Classical in
/-- ★★★★★★★★★★**引き戻した計量の特徴づけ**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

    `(f^* m).h V' (pullTrivOfBase f L W e_W) p = m.h W e_W (p ≫ f)`

★★これが `Classical.choice` の影を消す——`pullChartH_indep` により値はチャートに
依らないので、この等式は**すべての** `W`, `e_W`, `V' ≤ f⁻¹W` について成り立つ。
★★★`TorsorMetric.base` が落ちたのはまさにここであった。 -/
theorem pullback_h_pullTrivOfBase {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    (hL : IsLocallyTrivial Y L) (m : LocalMetric Y L) (W : Y.Opens)
    (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map f.base).obj W)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V' = ⊤) :
    (LocalMetric.pullback f hL m).h V' (pullTrivOfBase f L W eW hV'W) p
      = m.h W eW (p ≫ f) := by
  have hc : Nonempty (PullChart f L V' p) := nonempty_pullChart f hL V' p hp
  show (if hc : Nonempty (PullChart f L V' p) then
      pullChartH f m hc.some (pullTrivOfBase f L W eW hV'W) else 1) = _
  rw [dif_pos hc]
  have hself : pullChartH f m
      (⟨W, V', le_refl V', hV'W, hp, eW⟩ : PullChart f L V' p)
      (pullTrivOfBase f L W eW hV'W) = m.h W eW (p ≫ f) := by
    show ‖evalOn p V' hp (transUnit ((pullbackPre f).obj L) V'
        (pullTrivOfBase f L W eW hV'W)
        (trivialOfLe ((pullbackPre f).obj L) (le_refl V')
          (pullTrivOfBase f L W eW hV'W)))‖⁻¹ * _ = _
    rw [trivialOfLe_refl ((pullbackPre f).obj L) (pullTrivOfBase f L W eW hV'W),
      transUnit_self, evalOn_one, norm_one, inv_one, one_mul]
  rw [← hself]
  exact pullChartH_indep f m _ _ (pullTrivOfBase f L W eW hV'W)

/-! ## ★★★★★★★★★★算術直線束の引き戻し -/

/-- ★★★★★★★★★★**算術直線束の引き戻し** `φ^* L̄`。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★これが原文の「evident notion of pull-back」の中身である。
★★層は `pullbackPre`、局所自明性は `isLocallyTrivial_pullbackPre`、
計量は `LocalMetric.pullback`。 -/
noncomputable def AMetricPullback {X Y : Scheme.{0}} (f : X ⟶ Y) (L : AMetric Y) : AMetric X :=
  ⟨(pullbackPre f).obj L.sheaf, isLocallyTrivial_pullbackPre f L.sheaf L.triv,
    LocalMetric.pullback f L.triv L.metric⟩

/-! ### ★出典の紐付け(`.src`) -/

def PullChart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻しの計量のチャート)",
    sectionId := "genell-def-1-1-ii" }

def pullChartH_indep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(候補値がチャートに依らないこと)",
    sectionId := "genell-def-1-1-ii" }

def LocalMetric.pullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(計量の引き戻し)",
    sectionId := "genell-def-1-1-ii" }

def pullback_h_pullTrivOfBase.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻した計量の特徴づけ)",
    sectionId := "genell-def-1-1-ii" }

def AMetricPullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(算術直線束の引き戻し＝原文の evident notion of pull-back)",
    sectionId := "genell-def-1-1-ii" }

def AMetricPullback.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isLocallyTrivial_pullbackPre(局所自明性が引き戻しで保たれること)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isLocallyTrivial_pullbackPre") 3,
    .citation "[ABC3]" "transUnit_pullTrivOfBase(遷移単元は引き戻しと可換、§9-740)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.transUnit_pullTrivOfBase") 3,
    .citation "[ABC3]" "transUnit_pullTrivOfBase_shrinkW(W の縮小に依らないこと、§9-741)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.transUnit_pullTrivOfBase_shrinkW") 3,
    .implicitStep
      ("★原文は引き戻しを **evident** で畳んでいる。畳まれていたのは 6 段: " ++
       "層の引き戻し・局所自明性の保存・自明化の輸送・遷移単元の可換性(同じ W)・" ++
       "W の縮小に依らないこと・計量の構成") 3,
    .implicitStep
      ("★★残っている段: 引き戻しが**テンソル積と可換**であること" ++
       "(f^*(L̄ ⊗ M̄) ≅ f^*L̄ ⊗ f^*M̄ が等長)。" ++
       "これが ht_L(x) = deg_F(x_F^* L̄) の加法性に要る") 3 ]

end ABC3.Found.Arakelov
