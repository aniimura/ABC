import ABC3.Found.PGC.AxEpsilonDecay
import ABC3.Found.PGC.PGroupDescentToAxWild

/-!
# [pGC] `AxWildDescentDecay`（`θ ≤ 1` のつまみ）は**つまみではない** —— `AxLemma` と同値（循環）

## ★配られた問い「非空虚性」への答え

木は `Found/PGC/AxEpsilonDecay.lean:600-604` で
「★`C = p^{1/(p−1)}` を入れると Ax の定数 `p^{p/(p−1)²}` より真に良いので、
★`θ ≤ 1` は一般には成り立たないと見るべきである。★**非空虚性は測っていない。**」
と書いていた。★**測った。結果は木の見立てより強い。**

**(a) 警告の数の突き合わせは正しい**（§4 `axDecay_one_lt_axConstant`）。
`axDecay p 1 = p^{1/(p−1)} < p^{p/(p−1)²} = axConstant p`（`p ≥ 2`）。
指数を通分すると `(p−1)/(p−1)² < p/(p−1)²`、すなわち `p−1 < p` である。

**(b) ★★★しかし本当の理由は「良すぎる」ことではない —— `AxWildDescentDecay` は
`AxLemma` と同値である（循環）**（§2 `axWildDescentDecay_iff_axLemma`）。

`AxLemma K C` から `AxWildDescentDecay K C θ` が★**どんな `θ ≥ 0` でも**出る:
`x′ := algebraMap y`（`AxLemma` が返す `K` の元）を取ればよい。このとき

* `wildDepth K x′ = 0 < wildDepth K x`（`p ∣ deg minpoly x` だから）、
* `‖x − x′‖ ≤ C·ε`（`AxLemma` そのもの）、
* `‖σ • x′ − x′‖ = 0 ≤ θ·ε`（★`K` の元は `Γ_K` で固定される。`smul_algebraMap` 1 行）。

逆は木の `axLemma_of_wildDescentDecay`。⇒ `0 ≤ θ ≤ 1` の範囲で
★★**`AxWildDescentDecay K C θ ⟺ AxLemma K C`**、しかも
★**`θ` の値は何も効いていない**（§5 `axWildDescentDecay_inert`）。

⇒ ★★★**非空虚性の答え: `AxWildDescentDecay K C θ`（`0 ≤ θ ≤ 1`）が空虚でないのは
ちょうど `AxLemma K C` が真であるときである。掘っても何も減らない。**
`C = axDecay p 1` を入れると★**Ax の定理そのものの改良**を証明することになる
（§5 `axLemma_axConstant_of_axWildDescentDecay_axDecay_one`）。

★木の警告は結論として正しかったが、★理由は木が書いたもの（良すぎる）ではなく**循環**である。

## ★もう一方の道（`AxWildDescent K c`）は循環では**ない**

同じ議論は `AxWildDescent K c` にも半分だけ効く（§3 `axWildDescent_of_axLemma`）:
`AxLemma K C` は `∀ k, C ≤ c k` のとき `AxWildDescent K c` を与える。
★しかし `c = axDecay p` では `C > 1` に対して `∀ k, C ≤ axDecay p k` が**偽**なので
（前波 `PGroupDescentToAxWild.not_forall_le_axDecay`）、
★★**`AxWildDescent K (axDecay p)` は `AxLemma K C`（`C > 1`）からは出ない ⇒ 循環ではない**
（§6 `axDecay_route_not_circular`）。⇒ ★掘る価値があるのはこちらである。

## ★平均化の道（持ち場の証拠 3）—— 木がすでに測っていた

`NormalizedTraceDescent.lean` を読んだ（★本連鎖は本波で初めて開いた）:

* 跡の平均化は sharp な指数の★**ちょうど `(p−1)` 倍**を払うので `p ≥ 3` では閉じない
  （同 `:262 traceLoss_eq_sharpLoss_rpow`。`p = 2` なら一致し閉じる）。
* ★★同ファイル `:379 FirstJumpRoute.axLemma_of_firstJump` が★**閉じた十分条件**を持つ:
  `c k ≤ p^{j_k/(m_k·e_k)}`・`p^{k−1} ≤ m_k`・`(p−1)·j_k ≤ e_k` ⇒ `AxLemma K (axConstant p)`。
* ★この `c` は `k` に依って 1 に収束するので★**前波の②一様定数 no-go には当たらない**。
  ★また本波の前の波の③幾何減衰 no-go は `ε` の**変位の比**についての主張であり、
  こちらは**損失**についての主張なので★**当たらない**（別の量である）。
  ⇒ ★★`FirstJumpRoute` の仮説は 3 つの穴のどれにも塞がれていない。

★★**したがって「埋めなくてよい穴」が 1 つ見つかった**:
`AxWildDescentDecay` の道は掘る必要がない（循環だから）。
★残る 1 点は `FirstJumpRoute` が仮定として残している `AxWildDescent K c` を
「`c k ≤ p^{j_k/(m_k e_k)}`・`p^{k−1} ≤ m_k`」の形で埋めることである。
★ただし `AxWildDescent` は損失と `ε` の伸びを**同じ `c`** で押さえるので、
★★**`ε` の伸びを `c_k → 1` にできるか**が本当の 1 点である（★本波では測っていない）。

## ★測定器（§3 `le_of_axLemma`）

`d(x,K) ≥ D·ε` を実現する `x` が 1 つでもあれば `D ≤ C`。
⇒ ★`AxLemma K (axDecay p 1)` は**反例 1 つで死ぬ**。上の同値により、
その反例はそのまま `AxWildDescentDecay K (axDecay p 1) θ` の反例でもある。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は `AxEpsilonDecay.lean` と同じ項目（pGC 物理 p.6 Corollary 3.1）を指す。
   本ファイルの内容に対応する原典の文は無い。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★**抽象核に切り出さなかった。理由**: §1–§2 の中身は
   「`K` の元は `Γ_K` で固定される」という `K.closure` に固有の事実（`smul_algebraMap`）
   1 つだけで、切り出すと核が空になる。★切り出せなかった理由として記録する。
-/

namespace ABC3.Found.PGC

namespace AxKnob

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## §1 `K` の元は動かない -/

section Fixed

/-- ★`K` の元は `Γ_K` で固定されるので、変位は `0`。★これが §2 の循環の中身である。 -/
theorem norm_smul_algebraMap_sub_le (K : PAdicLocalField p) (y : K.carrier) {c : ℝ}
    (hc : 0 ≤ c) (σ : K.absGal) :
    ‖σ • algebraMap K.carrier K.closure y - algebraMap K.carrier K.closure y‖ ≤ c := by
  have h : σ • algebraMap K.carrier K.closure y = algebraMap K.carrier K.closure y :=
    smul_algebraMap σ y
  rw [h, sub_self, norm_zero]
  exact hc

/-- `p ∣ deg minpoly x` なら `0 < wildDepth K x`。 -/
theorem zero_lt_wildDepth_of_dvd (K : PAdicLocalField p) {x : K.closure}
    (hdvd : p ∣ (minpoly K.carrier x).natDegree) : 0 < wildDepth K x := by
  rcases Nat.eq_zero_or_pos (wildDepth K x) with h | h
  · exact absurd hdvd ((wildDepth_eq_zero_iff K x).mp h)
  · exact h

end Fixed

/-! ## §2 ★★★循環 —— `AxLemma` から `AxWildDescentDecay` が無料で出る -/

section Circular

/-- ★★★**循環の向き** —— `AxLemma K C` から `AxWildDescentDecay K C θ` が
★**どんな `θ ≥ 0` でも**出る（`x′` を `K` の中に取ればよい）。

★これで `AxEpsilonDecay.lean:605` の「`ε` が増えない形」は**つまみではなくなる**。 -/
theorem axWildDescentDecay_of_axLemma (K : PAdicLocalField p) {C θ : ℝ} (hθ : 0 ≤ θ)
    (h : AxLemma K C) : AxWildDescentDecay K C θ := by
  intro ε hε x hx hdvd
  obtain ⟨y, hy⟩ := h ε hε x hx
  refine ⟨algebraMap K.carrier K.closure y, ?_, hy, ?_⟩
  · rw [wildDepth_algebraMap]
    exact zero_lt_wildDepth_of_dvd K hdvd
  · exact fun σ => norm_smul_algebraMap_sub_le K y (mul_nonneg hθ hε) σ

/-- ★★★★**本ファイルの主結果** —— `1 ≤ C`・`0 ≤ θ ≤ 1` のとき
`AxWildDescentDecay K C θ ⟺ AxLemma K C`。

⇒ ★**非空虚性は `AxLemma K C` の真偽そのもの**である。掘っても何も減らない。 -/
theorem axWildDescentDecay_iff_axLemma (K : PAdicLocalField p) {C θ : ℝ}
    (hC : 1 ≤ C) (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) :
    AxWildDescentDecay K C θ ↔ AxLemma K C :=
  ⟨fun h => axLemma_of_wildDescentDecay K hC hθ0 hθ1 h,
   fun h => axWildDescentDecay_of_axLemma K hθ0 h⟩

/-- ★`AxLemma K C` は `∀ k, C ≤ c k` のとき `AxWildDescent K c` を与える。
★逆は出ない（積の条件が要る）ので、`AxWildDescent` は**循環ではない**。§6 を見よ。 -/
theorem axWildDescent_of_axLemma (K : PAdicLocalField p) {C : ℝ} {c : ℕ → ℝ}
    (hC0 : 0 ≤ C) (hc : ∀ k, C ≤ c k) (h : AxLemma K C) : AxWildDescent K c := by
  intro ε hε x hx hdvd
  obtain ⟨y, hy⟩ := h ε hε x hx
  refine ⟨algebraMap K.carrier K.closure y, ?_, ?_, ?_⟩
  · rw [wildDepth_algebraMap]
    exact zero_lt_wildDepth_of_dvd K hdvd
  · exact hy.trans (mul_le_mul_of_nonneg_right (hc _) hε)
  · exact fun σ =>
      norm_smul_algebraMap_sub_le K y (mul_nonneg (le_trans hC0 (hc _)) hε) σ

end Circular

/-! ## §3 `AxLemma` の単調性と測定器 -/

section Measure

/-- `AxLemma` は定数について単調。 -/
theorem axLemma_mono (K : PAdicLocalField p) {C C' : ℝ} (hCC : C ≤ C') (h : AxLemma K C) :
    AxLemma K C' := by
  intro e he x hx
  obtain ⟨y, hy⟩ := h e he x hx
  exact ⟨y, hy.trans (mul_le_mul_of_nonneg_right hCC he)⟩

/-- ★★**測定器** —— `d(x,K) ≥ D·ε` を実現する `x` が 1 つでもあれば `D ≤ C`。
⇒ ★`AxLemma K C` は反例 1 つで死ぬ。§2 の同値により
その反例は `AxWildDescentDecay K C θ` の反例でもある。 -/
theorem le_of_axLemma (K : PAdicLocalField p) {C : ℝ} (h : AxLemma K C) {e D : ℝ} (he : 0 < e)
    {x : K.closure} (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ e)
    (hlow : ∀ y : K.carrier, D * e ≤ ‖x - algebraMap K.carrier K.closure y‖) : D ≤ C := by
  obtain ⟨y, hy⟩ := h e (le_of_lt he) x hx
  have := (hlow y).trans hy
  exact le_of_mul_le_mul_right (by linarith) he

end Measure

/-! ## §4 木の警告の数の突き合わせ（正しい） -/

section Numeric

/-- ★木の警告（`AxEpsilonDecay.lean:600-603`）の数の突き合わせ:
`axDecay p 1 = p^{1/(p−1)} < p^{p/(p−1)²} = axConstant p`。★正しい。 -/
theorem axDecay_one_lt_axConstant (hp : 2 ≤ p) : axDecay p 1 < axConstant p := by
  have h2 : (2:ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
  have h1 : (1:ℝ) < (p : ℝ) := by linarith
  have hsq : (0:ℝ) < ((p:ℝ) - 1) ^ 2 := by nlinarith
  have hne : ((p:ℝ) - 1) ≠ 0 := ne_of_gt (by linarith)
  have key : (p:ℝ) / ((p:ℝ) - 1) ^ 2 - 1 / ((p:ℝ) - 1) = 1 / ((p:ℝ) - 1) ^ 2 := by
    field_simp
    ring
  have hpos : (0:ℝ) < 1 / ((p:ℝ) - 1) ^ 2 := by positivity
  rw [axDecay_one, axConstant]
  refine (Real.rpow_lt_rpow_left_iff h1).mpr ?_
  linarith

end Numeric

/-! ## §5 ★`θ` は効いていない -/

section Inert

/-- ★`0 ≤ θ ≤ 1` の範囲では `θ` の値は**何も効いていない**（全部同値）。 -/
theorem axWildDescentDecay_inert (K : PAdicLocalField p) {C θ θ' : ℝ} (hC : 1 ≤ C)
    (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) (hθ0' : 0 ≤ θ') (hθ1' : θ' ≤ 1) :
    AxWildDescentDecay K C θ ↔ AxWildDescentDecay K C θ' :=
  (axWildDescentDecay_iff_axLemma K hC hθ0 hθ1).trans
    (axWildDescentDecay_iff_axLemma K hC hθ0' hθ1').symm

theorem axWildDescentDecay_axDecay_one_iff (K : PAdicLocalField p) {θ : ℝ}
    (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) :
    AxWildDescentDecay K (axDecay p 1) θ ↔ AxLemma K (axDecay p 1) :=
  axWildDescentDecay_iff_axLemma K (one_le_axDecay p 1) hθ0 hθ1

/-- ★`C = axDecay p 1` で掘ると、出るのは `AxLemma K (axConstant p)` より**真に強い**
`AxLemma K (axDecay p 1)` である。★すなわち Ax の定理そのものの改良になる。 -/
theorem axLemma_axConstant_of_axWildDescentDecay_axDecay_one (K : PAdicLocalField p)
    (hp : 2 ≤ p) {θ : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1)
    (h : AxWildDescentDecay K (axDecay p 1) θ) : AxLemma K (axConstant p) :=
  axLemma_mono K (le_of_lt (axDecay_one_lt_axConstant hp))
    ((axWildDescentDecay_axDecay_one_iff K hθ0 hθ1).mp h)

end Inert

/-! ## §6 ★対比 —— `axDecay` の道は循環では**ない** -/

section Contrast

/-- ★★**対比** —— `AxWildDescentDecay`（`θ ≤ 1`）は `AxLemma` と同値（循環）だが、
`AxWildDescent K (axDecay p)` は `AxLemma K C`（`C > 1`）からは出ない。
⇒ ★**掘る価値があるのは後者である。** -/
theorem axDecay_route_not_circular (K : PAdicLocalField p) {C : ℝ} (hC : 1 < C) :
    (¬ (∀ k : ℕ, C ≤ axDecay p k))
      ∧ (∀ θ : ℝ, 0 ≤ θ → θ ≤ 1 → (AxWildDescentDecay K C θ ↔ AxLemma K C)) :=
  ⟨PGroupToAxWild.not_forall_le_axDecay hC,
   fun _ h0 h1 => axWildDescentDecay_iff_axLemma K (le_of_lt hC) h0 h1⟩

end Contrast

/-! ## §7 `.src` と 使っている公理の一覧 -/

section Src

def axWildDescentDecay_of_axLemma.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def axWildDescentDecay_iff_axLemma.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def le_of_axLemma.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def axDecay_one_lt_axConstant.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def axDecay_route_not_circular.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms norm_smul_algebraMap_sub_le
#print axioms zero_lt_wildDepth_of_dvd
#print axioms axWildDescentDecay_of_axLemma
#print axioms axWildDescentDecay_iff_axLemma
#print axioms axWildDescent_of_axLemma
#print axioms axLemma_mono
#print axioms le_of_axLemma
#print axioms axDecay_one_lt_axConstant
#print axioms axWildDescentDecay_inert
#print axioms axWildDescentDecay_axDecay_one_iff
#print axioms axLemma_axConstant_of_axWildDescentDecay_axDecay_one
#print axioms axDecay_route_not_circular

end AxKnob

end ABC3.Found.PGC
