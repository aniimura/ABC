import ABC3.Found.PGC.NormalizedTraceDescent
import ABC3.Found.PGC.WildDescentDistanceOnly

/-!
# [pGC] `FirstJumpRoute` は `axLemma_of_axDecay` の**言い換え**だった —— 新しい道ではない

## ★配られた問いへの答え

**(1) `hjump` の係数（`e` か `p^{k+1}·e` か）—— ★別のものである。同じ定理を別の体に当てている。**

`NormalizedTraceDescent.lean:214-221` の docstring が定義を書いている:
`j` = ★**下の層** `E₁/F` の跳び、`m = e(L/E₁)`、`e = e_{E₁}`、`m·e = e_L`。
⇒ `hjump : (p−1)·j ≤ e` は `RamificationJumpBound.lean:334
sub_one_mul_le_of_norm_natCast_eq_pow`（`(p−1)·i ≤ e_L`、`e_L = v_L(p)`）を
★**下の層 `E₁`** に当てたもの。
一方 `harith` (4)（`WildBreakPowDegree.hupper_of_totallyRamified_pow`）は
`(p−1)·u_k ≤ p^{k+1}·e` で、★**塔の頂上 `M`** に当てたもの（`p^{k+1}e = v_M(p) = e_M`）。
⇒ ★**同じ古典的評価を違う体に当てた 2 つの主張**であって、矛盾していない。
★原典で「どちらが正しいか」を測る必要は無い（どちらも同じ 1 本の帰結である）。

**(2) ★★★しかし `FirstJumpRoute` の 4 つの仮説は `c k ≤ axDecay p k` と**同値**である**
（§1 `exists_firstJump_data_iff`）。

* `⇒` は木の `JumpArith.rpow_div_le_axDecay` そのもの。
* `⇐` は `j = 1`・`m = p^{k−1}`・`e = p−1` を取れば★**等号**で実現する
  （`axDecay_eq_rpow_one_div`）。

⇒ ★★`FirstJumpRoute.axLemma_of_firstJump` は `axLemma_of_axDecay` の**再パラメータ化**で、
★**新しい情報を 1 ビットも持っていない**。
木の docstring「★★これが `AxWildDescent K (axDecay p)` を出す**閉じた十分条件**である」は
字面としては真だが、★同値なので「閉じた」という語が示唆するほどの前進ではない。

**(3) ★★具体層は木がすでに測っていて、破れている。**

`WildDescentDistanceOnly.lean`（★本連鎖は本波で初めて開いた）は
`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`（`k = 2`, `e_L = 18`）で★**実在する `x`** について `i − gain = 4` を
厳密整数演算で出しており、要求 `(p−1)·p^{k−1}·(i−gain) ≤ e_L` すなわち
`2·3·4 = 24 ≤ 18` は**偽**である。
本ファイルは同じ破れを `FirstJumpRoute` の仮説の言葉に書き直した（§2）:

★`not_exists_firstJump_data_cyclotomic` —— `3 ≤ m`・`8 ≤ e`・`m·e = 18` を満たす
`(m, e)` は**存在しない**（`3·8 = 24 > 18`）。
⇒ ★★**`FirstJumpRoute` の `hform`/`hm`/`hjump` はこの測定点では満たせない。**

## ★★★自分の前波の主張の訂正（2 つとも外れていた）

1. 前波で「残る 1 点は★**`ε` の伸びを `c_k → 1` にできるか**」と書いたが★**外れ**。
   `WildDescentDistanceOnly.lean:201 axWildDescent_of_dist` が
   ★**`ε` の伸びは `1 ≤ c k` と超距離性から無償**であることを既に定理にしている
   （`σx′−x′ = σ(x′−x) + (σx−x) + (x−x′)` を `max` で潰すだけ）。
   ⇒ ★残っているのは**距離の評価 1 本**だけ。
   §3 でそれを使い★**「距離の評価 1 本から `AxSenTate K`」**を 1 本の定理にした
   （`axSenTate_of_dist_axDecay` / `axLemma_of_dist_axDecay`）。
2. 前波で「`FirstJumpRoute` の仮説は 3 つの穴のどれにも塞がれていない」と書いたが
   ★**測り方が甘かった**。②一様定数 no-go には確かに当たらないが、
   `WildDescentDistanceOnly` の測定 3 が★**同じ不等式**を実データで破っている。
   ⇒ ★★これが**4 つ目の穴**である（★ただし木自身が「反例ではなく**反証候補**」と
   明記している —— `M` の外の `x′` を排除していない）。

## ★木の 2 つの docstring が食い違っている（名指しで記録。どちらも書き換えない）

* `NormalizedTraceDescent.lean:374-386`:
  「★★**定数の勘定はここで完全に閉じている**。残っているのは『そういう `x′` を作る』具体層だけ」
* `WildDescentDistanceOnly.lean:10-16`:
  「★★**`AxWildDescent K (axDecay p)` は閉じていない。**」

★**後者が新しく、かつ厳密整数演算の裏づけを持つ。**
本ファイルの `not_exists_firstJump_data_cyclotomic` は、後者の数値が前者の仮説を
**満たさない**ことを `ℕ` の 1 行で示す。
⇒ ★前者の「定数の勘定は閉じている」は**同値だから**真だが、
★**「あとは具体層だけ」という見立ては後者に否定されている**。

## ★①不分岐側との関係

★本波も独立である。上の破れは `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`（完全分岐）で起きており、
`TotallyRamifiedLayer.lean:342 not_valK_of_norm_eq`（不分岐側）とは別の穴である。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は上流と同じ項目（pGC 物理 p.6 Corollary 3.1）。本ファイルの内容に
   対応する原典の文は無い。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. §1 は `ℝ` と `ℕ` だけ、§2 は `ℕ` だけで、分岐・付値・Galois の語彙が 1 語も出ない。
   §3 だけが `PAdicLocalField` に触る。
-/

namespace ABC3.Found.PGC

namespace FirstJumpEquiv

/-! ## §1 ★★★`FirstJumpRoute` の仮説は `c k ≤ axDecay p k` と同値 -/

section Core

variable {p : ℕ} [Fact p.Prime]

/-- `axDecay p k = p^{1/(p^{k−1}·(p−1))}`（`ℕ` の指数で書いた形）。
★§1 の逆向き（`j = 1`, `m = p^{k−1}`, `e = p−1` で等号）に使う。 -/
theorem axDecay_eq_rpow_one_div (k : ℕ) :
    axDecay p k
      = (p : ℝ) ^ ((1 : ℝ) / (((p ^ (k - 1) : ℕ) : ℝ) * (((p - 1 : ℕ)) : ℝ))) := by
  have hp1 : (1:ℕ) ≤ p := (Fact.out : p.Prime).one_lt.le
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  rw [axDecay]
  congr 1
  rw [Nat.cast_pow, Nat.cast_sub hp1, Nat.cast_one, one_div_pow, one_div, one_div,
    ← mul_inv, inv_eq_one_div, mul_comm]

/-- ★★★★**本ファイルの主結果** —— `NormalizedTraceDescent.lean:379
FirstJumpRoute.axLemma_of_firstJump` の 4 つの仮説（`hform`/`hm`/`he`/`hjump`）は
★**`c k ≤ axDecay p k` と同値**である。

⇒ ★`FirstJumpRoute` は `axLemma_of_axDecay` の**再パラメータ化**であって、
★**新しい情報を 1 ビットも持っていない**。 -/
theorem exists_firstJump_data_iff {k : ℕ} {ck : ℝ} :
    (∃ j m e : ℕ, ck ≤ (p : ℝ) ^ ((j : ℝ) / ((m : ℝ) * (e : ℝ))) ∧
        p ^ (k - 1) ≤ m ∧ 0 < e ∧ (p - 1) * j ≤ e)
      ↔ ck ≤ axDecay p k := by
  constructor
  · rintro ⟨j, m, e, hform, hm, he, hjump⟩
    exact hform.trans (JumpArith.rpow_div_le_axDecay hm he hjump)
  · intro h
    have hp2 : 2 ≤ p := (Fact.out : p.Prime).two_le
    refine ⟨1, p ^ (k - 1), p - 1, ?_, le_rfl, by omega, by omega⟩
    rw [Nat.cast_one, ← axDecay_eq_rpow_one_div (p := p) k]
    exact h

end Core

/-! ## §2 ★★測定点（`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`）では仮説が満たせない -/

section Concrete

/-- ★★**測定点では仮説が満たせない** —— `WildDescentDistanceOnly.lean` の
`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`（`k = 2`, `e_L = 18`, 実在する `x` で `j = i − gain = 4`）では
`3 ≤ m`・`8 ≤ e`・`m·e = 18` を満たす `(m, e)` が**存在しない**（`3·8 = 24 > 18`）。

★★これは `FirstJumpRoute` の docstring「あとは具体層だけ」への反証候補である
（★木自身が『反例ではなく反証候補』と明記している。`M` の外の `x′` は排除していない）。 -/
theorem not_exists_firstJump_data_cyclotomic :
    ¬ ∃ m e : ℕ, 3 ^ (2 - 1) ≤ m ∧ 0 < e ∧ (3 - 1) * 4 ≤ e ∧ m * e = 18 := by
  rintro ⟨m, e, hm, _he, hjump, hme⟩
  have h1 : 3 ≤ m := by simpa using hm
  have h2 : 8 ≤ e := by omega
  have : 3 * 8 ≤ m * e := Nat.mul_le_mul h1 h2
  omega

end Concrete

/-! ## §3 ★残っているのは「距離の評価 1 本」だけ -/

section Route

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- §1 を全ての `k` について束ねた形。 -/
theorem forall_firstJump_data_iff {c : ℕ → ℝ} :
    (∀ k : ℕ, ∃ j m e : ℕ, c k ≤ (p : ℝ) ^ ((j : ℝ) / ((m : ℝ) * (e : ℝ))) ∧
        p ^ (k - 1) ≤ m ∧ 0 < e ∧ (p - 1) * j ≤ e)
      ↔ (∀ k : ℕ, c k ≤ axDecay p k) :=
  forall_congr' fun _ => exists_firstJump_data_iff

/-- ★★★**鎖の現在地を 1 本にまとめた形** —— ★**距離の評価 1 本から `AxSenTate K`**。

`ε` の伸びは `WildDescentDistanceOnly.lean:201 axWildDescent_of_dist` が無償で落とす
（★前波に本連鎖が「残る 1 点は `ε` の伸び」と書いたのは**外れ**だった）。
⇒ ★残っているのは「wild 深さ `k` の `x` に対し `‖x − x′‖ ≤ axDecay p k · ε` を満たす
深さ `< k` の `x′` を作る」ただ 1 本である。 -/
theorem axSenTate_of_dist_axDecay (K : PAdicLocalField p)
    (h : ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
      p ∣ (minpoly K.carrier x).natDegree →
        ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧
          ‖x - x'‖ ≤ axDecay p (wildDepth K x) * ε) :
    AxSenTate K :=
  axSenTate_of_axDecay K (one_le_axDecay p) (fun _ => le_rfl)
    (axWildDescent_of_dist K (one_le_axDecay p) h)

/-- 同じものの `AxLemma K (axConstant p)` 版。 -/
theorem axLemma_of_dist_axDecay (K : PAdicLocalField p)
    (h : ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
      p ∣ (minpoly K.carrier x).natDegree →
        ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧
          ‖x - x'‖ ≤ axDecay p (wildDepth K x) * ε) :
    AxLemma K (axConstant p) :=
  axLemma_of_axDecay K (one_le_axDecay p) (fun _ => le_rfl)
    (axWildDescent_of_dist K (one_le_axDecay p) h)

end Route

/-! ## §4 `.src` と 使っている公理の一覧 -/

section Src

def exists_firstJump_data_iff.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def not_exists_firstJump_data_cyclotomic.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def axSenTate_of_dist_axDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms axDecay_eq_rpow_one_div
#print axioms exists_firstJump_data_iff
#print axioms not_exists_firstJump_data_cyclotomic
#print axioms forall_firstJump_data_iff
#print axioms axSenTate_of_dist_axDecay
#print axioms axLemma_of_dist_axDecay

end FirstJumpEquiv

end ABC3.Found.PGC
