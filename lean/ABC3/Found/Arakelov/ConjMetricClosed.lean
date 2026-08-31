/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.ConjMetric
import ABC3.Found.Arakelov.TensorMetric
import ABC3.Found.Arakelov.PullbackMetric
import ABC3.Found.Arakelov.AMetricIso
import ABC3.Found.Arakelov.AMetricPic
import ABC3.Found.Arakelov.PullbackPic
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★共役両立性は `1`・`⊗`・引き戻し・等長で**閉じている**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> that is compatible [in the evident sense] with the complex conjugation automorphism ιX.

## ★★★★★★★★★★これは何か——`APic(X)` が**群になる**ということ

原文は Definition 1.1 (i) で

> The isomorphism classes of arithmetic line bundles on X, together with the operation
> of tensor product, thus determine a group APic(X).

と書く。★ここでの「arithmetic line bundle」は**共役両立な計量を持つもの**である。
★★だから `APic(X)` が群になるには、共役両立性が

| 演算 | 本ファイルの定理 |
|---|---|
| 単位元 `Ō_X` | `AMetric.one_isConjCompatible` |
| テンソル積 `⊗` | `AMetric.mul_isConjCompatible` |
| 逆元 | `AInv.inv_isConjCompatible` |
| 等長同型（＝類の代表の取り替え） | `AMetric.isConjCompatible_of_isometric` |
| 引き戻し `φ^*`（(ii) 用） | `AMetricPullback_isConjCompatible` |

のすべてで閉じている必要がある。★★★本ファイルはそれを示し、
**`APicC X : Subgroup (APicM X)`** として結晶させる。

## ★★★★★★★逆元が自動で共役両立になること

★これは自明ではない。`AInv` は逆元を**データとして持つ**ので、
「与えられた逆元がたまたま共役両立でない」可能性を潰す必要がある。

    `L̄ ⊗ L̄⁻¹ ≅ Ō_X`（等長）  ⇒  `L̄ ⊗ L̄⁻¹` は共役両立
    `h_{A⊗B} = h_A · h_B`（`IsTensorOf`）  ⇒  `h_A > 0` で割れる

★★これが `LocalMetric.isConjCompatible_of_tensor` である
——**因子の一方が共役両立なら、もう一方も共役両立**。

## ★★★★★機構——チャートを `ι_X p` へ運ぶ

`LocalMetric.tensor` と `LocalMetric.pullback` の `h` は
`Nonempty (…Chart …)` の**選択**で書かれているので、
`p` の側のチャートと `ι_X p` の側のチャートが**別物**になる。

★そこで `TensorChart.conj` / `PullChart.conj` でチャートを運び、
在庫の `chartH_indep` / `pullChartH_indep`（値はチャートに依らない）で
選択の食い違いを吸収する。★★これが本ファイルの配管である。

## ★測定の記録

`pullChartH` は `m.h c.W c.eW (p ≫ f)` を含む。
★`conjPoint_comp`（`ι_X (p ≫ f) = (ι_X p) ≫ f`、`ArcConjInvol.lean`）で
底の側の共役に化けるので、底の計量の共役両立性がそのまま効く。
★★そのとき `(p ≫ f) ⁻¹ᵁ c.W = ⊤` が要る——`c.hV'W` と `c.hpV'` から出る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{0}}

/-! ## ★★★(1) 単位元 `Ō_X` -/

open scoped Classical in
/-- ★★★★**構造層の標準計量は `ι_X` と両立する**。

原文 (GenEll p.3):
> that is compatible [in the evident sense] with the complex conjugation automorphism ιX.

★`h` は遷移単元のノルムの逆数であり、`evalOn_conjPoint` で共役になるだけだからである。 -/
theorem structLocalMetric_isConjCompatible (X : Scheme.{0}) :
    (structLocalMetric X).IsConjCompatible := by
  intro V e p hp
  show (if h : conjPoint p ⁻¹ᵁ V = ⊤ then
      ‖evalOn (conjPoint p) V h (transUnit (𝟙_ X.PresheafOfModules) V (baseTriv X V) e)‖⁻¹
    else 1) = (if h : p ⁻¹ᵁ V = ⊤ then
      ‖evalOn p V h (transUnit (𝟙_ X.PresheafOfModules) V (baseTriv X V) e)‖⁻¹ else 1)
  rw [dif_pos (conjPoint_preimage_eq_top hp), dif_pos hp]
  congr 1
  rw [evalOn_conjPoint V p hp]
  exact RCLike.norm_conj _

/-! ## ★★★★★(2) テンソル積 -/

/-- ★チャートを `ι_X p` へ運ぶ（`p` と `ι_X p` は同じ開集合に落ちる）。 -/
def TensorChart.conj {A B : X.PresheafOfModules} {V : X.Opens}
    {p : Spec (CommRingCat.of ℂ) ⟶ X} (c : TensorChart A B V p) :
    TensorChart A B V (conjPoint p) :=
  ⟨c.W, c.hWV, conjPoint_preimage_eq_top c.hpW, c.eA, c.eB⟩

/-- ★★★**運んだチャートの上では候補値が変わらない**。 -/
theorem chartH_conj {A B : X.PresheafOfModules} {mA : LocalMetric X A} {mB : LocalMetric X B}
    (hA : mA.IsConjCompatible) (hB : mB.IsConjCompatible) {V : X.Opens}
    {p : Spec (CommRingCat.of ℂ) ⟶ X} (c : TensorChart A B V p)
    (f : (restrictPresheafFunctor X V).obj (A ⊗ B) ≅ 𝟙_ (PresheafModulesOn X V)) :
    chartH mA mB c.conj f = chartH mA mB c f := by
  show ‖evalOn (conjPoint p) c.W (conjPoint_preimage_eq_top c.hpW)
      (transUnit (A ⊗ B) c.W (tensorTriv c.eA c.eB) (trivialOfLe (A ⊗ B) c.hWV f))‖⁻¹
      * (mA.h c.W c.eA (conjPoint p) * mB.h c.W c.eB (conjPoint p)) = _
  rw [hA c.W c.eA p c.hpW, hB c.W c.eB p c.hpW]
  congr 2
  rw [evalOn_conjPoint c.W p c.hpW]
  exact RCLike.norm_conj _

open scoped Classical in
/-- ★★★★★★★★**テンソル積は共役両立性を保つ**。

原文 (GenEll p.3):
> that is compatible [in the evident sense] with the complex conjugation automorphism ιX.

★選択されたチャートの食い違いは `chartH_indep` で吸収する。 -/
theorem LocalMetric.tensor_isConjCompatible {A B : X.PresheafOfModules}
    (hlA : IsLocallyTrivial X A) (hlB : IsLocallyTrivial X B)
    {mA : LocalMetric X A} {mB : LocalMetric X B}
    (hA : mA.IsConjCompatible) (hB : mB.IsConjCompatible) :
    (LocalMetric.tensor hlA hlB mA mB).IsConjCompatible := by
  intro V f p hp
  have hc : Nonempty (TensorChart A B V p) := nonempty_tensorChart hlA hlB V p hp
  have hc' : Nonempty (TensorChart A B V (conjPoint p)) := ⟨hc.some.conj⟩
  show (if h : Nonempty (TensorChart A B V (conjPoint p)) then chartH mA mB h.some f else 1)
    = (if h : Nonempty (TensorChart A B V p) then chartH mA mB h.some f else 1)
  rw [dif_pos hc', dif_pos hc, chartH_indep mA mB hc'.some hc.some.conj f]
  exact chartH_conj hA hB hc.some f

/-! ## ★★★★★★(3) 引き戻し -/

/-- ★チャートを `ι_X p` へ運ぶ。 -/
def PullChart.conj {X Y : Scheme.{0}} {f : X ⟶ Y} {L : Y.PresheafOfModules}
    {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X} (c : PullChart f L V p) :
    PullChart f L V (conjPoint p) :=
  ⟨c.W, c.V', c.hV'V, c.hV'W, conjPoint_preimage_eq_top c.hpV', c.eW⟩

/-- ★★底の側の点 `p ≫ f` も `W` に落ちる。 -/
theorem PullChart.comp_preimage_eq_top {X Y : Scheme.{0}} {f : X ⟶ Y} {L : Y.PresheafOfModules}
    {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X} (c : PullChart f L V p) :
    (p ≫ f) ⁻¹ᵁ c.W = ⊤ := by
  rw [Scheme.Hom.comp_preimage]
  exact preimage_eq_top_of_le c.hV'W c.hpV'

/-- ★★★**運んだチャートの上では候補値が変わらない**。

★`m.h c.W c.eW (ι_X p ≫ f) = m.h c.W c.eW (ι_X (p ≫ f))`（`conjPoint_comp`）
なので、底の計量の共役両立性がそのまま効く。 -/
theorem pullChartH_conj {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    {m : LocalMetric Y L} (hm : m.IsConjCompatible) {V : X.Opens}
    {p : Spec (CommRingCat.of ℂ) ⟶ X} (c : PullChart f L V p)
    (e : (restrictPresheafFunctor X V).obj ((pullbackPre f).obj L)
      ≅ 𝟙_ (PresheafModulesOn X V)) :
    pullChartH f m c.conj e = pullChartH f m c e := by
  show ‖evalOn (conjPoint p) c.V' (conjPoint_preimage_eq_top c.hpV')
      (transUnit ((pullbackPre f).obj L) c.V' (pullTrivOfBase f L c.W c.eW c.hV'W)
        (trivialOfLe ((pullbackPre f).obj L) c.hV'V e))‖⁻¹
      * m.h c.W c.eW (conjPoint p ≫ f) = _
  rw [← conjPoint_comp f p, hm c.W c.eW (p ≫ f) c.comp_preimage_eq_top]
  congr 2
  rw [evalOn_conjPoint c.V' p c.hpV']
  exact RCLike.norm_conj _

open scoped Classical in
/-- ★★★★★★★★**引き戻しは共役両立性を保つ**。

原文 (GenEll p.3):
> that is compatible [in the evident sense] with the complex conjugation automorphism ιX. -/
theorem LocalMetric.pullback_isConjCompatible {X Y : Scheme.{0}} (f : X ⟶ Y)
    {L : Y.PresheafOfModules} (hL : IsLocallyTrivial Y L) {m : LocalMetric Y L}
    (hm : m.IsConjCompatible) :
    (LocalMetric.pullback f hL m).IsConjCompatible := by
  intro V e p hp
  have hc : Nonempty (PullChart f L V p) := nonempty_pullChart f hL V p hp
  have hc' : Nonempty (PullChart f L V (conjPoint p)) := ⟨hc.some.conj⟩
  show (if h : Nonempty (PullChart f L V (conjPoint p)) then pullChartH f m h.some e else 1)
    = (if h : Nonempty (PullChart f L V p) then pullChartH f m h.some e else 1)
  rw [dif_pos hc', dif_pos hc, pullChartH_indep f m hc'.some hc.some.conj e]
  exact pullChartH_conj f hm hc.some e

/-! ## ★★★★★★★(4) 等長同型で移ること -/

/-- ★★★★★★**共役両立性は等長同型で移る**。

★★これがあって初めて条件が `APic(X)` の**類の上**で意味を持つ。 -/
theorem AMetric.isConjCompatible_of_isIsometry {L M : AMetric X} {φ : L.sheaf ≅ M.sheaf}
    (hφ : IsIsometry L M φ) (hL : L.IsConjCompatible) : M.IsConjCompatible := by
  intro V e p hp
  rw [← hφ V e (conjPoint p) (conjPoint_preimage_eq_top hp), ← hφ V e p hp]
  exact hL V (pullTriv φ V e) p hp

theorem AMetric.isConjCompatible_of_isometric {L M : AMetric X} (h : Isometric L M)
    (hL : L.IsConjCompatible) : M.IsConjCompatible :=
  h.elim fun _ hφ => AMetric.isConjCompatible_of_isIsometry hφ hL

/-! ## ★★★★★★★★(5) 因子の一方から他方へ -/

/-- ★★★★★★★★**テンソル積とその一方の因子が共役両立なら、他方も共役両立**。

★機構は `IsTensorOf`（`h_{A⊗B} = h_A · h_B`）で割り算するだけ
——`h_A > 0` なので `mul_left_cancel₀` が効く。
★★両方が自明になる小さい `W` へ降り、`restrict` の欄で戻る。 -/
theorem LocalMetric.isConjCompatible_of_tensor {A B : X.PresheafOfModules}
    (hlA : IsLocallyTrivial X A) (hlB : IsLocallyTrivial X B)
    {mA : LocalMetric X A} {mB : LocalMetric X B} {m : LocalMetric X (A ⊗ B)}
    (ht : IsTensorOf mA mB m) (hm : m.IsConjCompatible) (hA : mA.IsConjCompatible) :
    mB.IsConjCompatible := by
  intro V eB p hp
  obtain ⟨c⟩ := nonempty_tensorChart hlA hlB V p hp
  set eB' := trivialOfLe B c.hWV eB with heB'
  have h1 : m.h c.W (tensorTriv c.eA eB') (conjPoint p)
      = mA.h c.W c.eA (conjPoint p) * mB.h c.W eB' (conjPoint p) :=
    ht c.W c.eA eB' (conjPoint p) (conjPoint_preimage_eq_top c.hpW)
  have h2 : m.h c.W (tensorTriv c.eA eB') p = mA.h c.W c.eA p * mB.h c.W eB' p :=
    ht c.W c.eA eB' p c.hpW
  have key : mA.h c.W c.eA p * mB.h c.W eB' (conjPoint p)
      = mA.h c.W c.eA p * mB.h c.W eB' p :=
    calc mA.h c.W c.eA p * mB.h c.W eB' (conjPoint p)
        = mA.h c.W c.eA (conjPoint p) * mB.h c.W eB' (conjPoint p) := by
          rw [hA c.W c.eA p c.hpW]
      _ = m.h c.W (tensorTriv c.eA eB') (conjPoint p) := h1.symm
      _ = m.h c.W (tensorTriv c.eA eB') p := hm c.W (tensorTriv c.eA eB') p c.hpW
      _ = mA.h c.W c.eA p * mB.h c.W eB' p := h2
  have hcancel := mul_left_cancel₀ (mA.pos c.W c.eA p).ne' key
  rw [← mB.restrict c.hWV eB (conjPoint p) (conjPoint_preimage_eq_top c.hpW),
    ← mB.restrict c.hWV eB p c.hpW]
  exact hcancel

/-! ## ★★★★★★★★★(6) 算術直線束の側 -/

/-- ★★★**単位元 `Ō_X` は共役両立**。 -/
theorem AMetric.one_isConjCompatible (X : Scheme.{0}) : (1 : AMetric X).IsConjCompatible :=
  structLocalMetric_isConjCompatible X

/-- ★★★★★★**テンソル積は共役両立性を保つ**。 -/
theorem AMetric.mul_isConjCompatible {L M : AMetric X} (hL : L.IsConjCompatible)
    (hM : M.IsConjCompatible) : (L * M).IsConjCompatible :=
  LocalMetric.tensor_isConjCompatible L.triv M.triv hL hM

/-- ★★★★★★**引き戻しは共役両立性を保つ**（Definition 1.1 (ii) 用）。 -/
theorem AMetricPullback_isConjCompatible {X Y : Scheme.{0}} (f : X ⟶ Y) {L : AMetric Y}
    (hL : L.IsConjCompatible) : (AMetricPullback f L).IsConjCompatible :=
  LocalMetric.pullback_isConjCompatible f L.triv hL

/-- ★★★★★★★**積と一方の因子から他方の因子へ**。 -/
theorem AMetric.isConjCompatible_of_mul {L M : AMetric X}
    (h : (L * M).IsConjCompatible) (hL : L.IsConjCompatible) : M.IsConjCompatible :=
  LocalMetric.isConjCompatible_of_tensor L.triv M.triv
    (isTensorOf_tensor L.triv M.triv L.metric M.metric) h hL

/-! ## ★★★★★★★★★★(7) `APic(X)` の中で群になること -/

/-- ★★★可逆な算術直線束が共役両立であること。 -/
def AInv.IsConjCompatible (L : AInv X) : Prop := L.carrier.IsConjCompatible

/-- ★★★★★★★★**逆元は自動で共役両立になる**。

★`L̄ ⊗ L̄⁻¹ ≅ Ō_X` が等長なので `L̄ ⊗ L̄⁻¹` は共役両立、
そこから `AMetric.isConjCompatible_of_mul` で `L̄⁻¹` が出る。
★★これが無いと「共役両立な類」は群にならない。 -/
theorem AInv.inv_isConjCompatible {L : AInv X} (h : L.IsConjCompatible) :
    L.inv.IsConjCompatible :=
  AMetric.isConjCompatible_of_mul
    (AMetric.isConjCompatible_of_isometric (isometric_symm L.isInv)
      (AMetric.one_isConjCompatible X)) h

/-- ★★★★**共役両立性は類の代表の取り替えで変わらない**。 -/
theorem AInv.isConjCompatible_congr {L M : AInv X} (h : (AInv.setoid X).r L M) :
    L.IsConjCompatible = M.IsConjCompatible :=
  propext ⟨fun hL => AMetric.isConjCompatible_of_isometric h hL,
    fun hM => AMetric.isConjCompatible_of_isometric (isometric_symm h) hM⟩

/-- ★★★★★★**共役両立性は同型類の上の性質である**。 -/
def APicM.IsConjCompatible : APicM X → Prop :=
  Quotient.lift AInv.IsConjCompatible (fun _ _ h => AInv.isConjCompatible_congr h)

/-- ★★★★★★★★★★**[GenEll] Definition 1.1, (i)** ——
**共役両立な算術直線束の同型類は `APic(X)` の部分群をなす**。

原文 (GenEll p.3):
> that is compatible [in the evident sense] with the complex conjugation automorphism ιX.

★★原文の `APic(X)` は「(共役両立な)算術直線束の同型類がテンソル積でなす群」である。
★★★本定義がその**群であること**を保証する
——単位元・積・逆元のすべてで共役両立性が閉じている。 -/
noncomputable def APicC (X : Scheme.{0}) : Subgroup (APicM X) where
  carrier := {a | APicM.IsConjCompatible a}
  mul_mem' := by
    rintro a b ha hb
    induction a using Quotient.ind with
    | _ L =>
      induction b using Quotient.ind with
      | _ M => exact AMetric.mul_isConjCompatible ha hb
  one_mem' := AMetric.one_isConjCompatible X
  inv_mem' := by
    rintro a ha
    induction a using Quotient.ind with
    | _ L => exact AInv.inv_isConjCompatible ha

/-- ★★★★★**引き戻しは可逆な算術直線束の共役両立性を保つ**。 -/
theorem AInv.pullback_isConjCompatible {X Y : Scheme.{0}} (f : X ⟶ Y) {L : AInv Y}
    (h : L.IsConjCompatible) : (AInv.pullback f L).IsConjCompatible :=
  AMetricPullback_isConjCompatible f h

/-! ## ★出典の紐付け(`.src`) -/

def LocalMetric.tensor_isConjCompatible.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(テンソル積は共役両立性を保つ)",
    sectionId := "genell-def-1-1-i" }

def LocalMetric.pullback_isConjCompatible.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻しは共役両立性を保つ)",
    sectionId := "genell-def-1-1-ii" }

def AMetric.isConjCompatible_of_isIsometry.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(共役両立性は等長同型で移る)",
    sectionId := "genell-def-1-1-i" }

def AInv.inv_isConjCompatible.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(逆元は自動で共役両立になる)",
    sectionId := "genell-def-1-1-i" }

def APicC.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(共役両立な算術直線束の同型類は APic(X) の部分群をなす)",
    sectionId := "genell-def-1-1-i" }

def APicC.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isTensorOf_tensor(h_{A⊗B} = h_A · h_B)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isTensorOf_tensor") 3,
    .citation "[ABC3]" "chartH_indep / pullChartH_indep(値はチャートに依らない)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.chartH_indep") 3,
    .citation "[ABC3]" "evalOn_conjPoint(ι_X は値を複素共役にする、§9-800)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.evalOn_conjPoint") 3,
    .implicitStep
      ("★原文は「The isomorphism classes of arithmetic line bundles on X, together with " ++
       "the operation of tensor product, thus determine a group APic(X).」と書くが、" ++
       "★★そこで畳まれていたのは「共役両立性が 1・⊗・逆元・等長で閉じていること」である") 3 ]

end ABC3.Found.Arakelov
