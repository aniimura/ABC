import ABC3.Found.PGC.DepthApproxToAxLemma

/-!
# [pGC] ★★★訂正: `Gained*` は「一気に」ではなく **1 層**下る ＋ `AxWildDescent` の第 3 条件は**無料**

## ★★★訂正（自己訂正 19 度目）—— 前波の私の読みが誤り

前波（第 1131）の私は本文にこう書いた（逐語）:

> `F = K.carrier` と取り、`‖τ x − x‖ ≤ ε` を使えば … これは **`AxLemma K (axConstant p)` そのもの**である。
> ⇒ ★★★**`AxWildDescent` も `axLemma_of_wildDescent_Icc` も基底段の別扱いも要らない。**

★★**誤り。** `Gained*` 族の `F` は ★**`[M : F] = p`**（`M` の**すぐ下**の体）である。
`F = K.carrier` と取れるのは `[M:K] = p`、すなわち ★**wild 深さ 1 のときだけ**。

証拠 2 つ（本波で読んだ）:

1. `GainedTowerStep.lean:404` の仮説（逐語）
   `(hdeg : (minpoly F π).natDegree = p)` と
   `(htop : Algebra.adjoin F ({π} : Set M) = ⊤)` ⇒ `[M : F] = p`。
2. `GainedTowerModel.lean:347` の docstring（逐語）
   「さらに ★**上の体 `F` 自身**(`= twr g p k`、**`[M:F] = p`**)と `hdeg` / `htop`。」
   結論も `∃ y : twr g p k` であり、`[M : twr g p k] = p`。

⇒ ★★`Gained*` 族は ★**1 層下る**（`M` → `F`）。★「一気に降りる」ではない。
★`PureStepSetup`（`k = 0` 専用）と同じ性質で、違うのは**払う値段**だけである。

## ★負の結果 —— 1 層あたり `∏_{j=1}^{k+1} axDecay p j` を払うと閉じない（§2、証明つき）

`Gained*` が 1 層に払うのは `∏_{j ∈ Icc 1 (k+1)} axDecay p j`（`k+1` = wild 深さ）。
★これを深さごとの `c k` として `AxWildDescent` に積むと、
`∏_{k ∈ Icc 1 n} (∏_{j ∈ Icc 1 k} axDecay p j)` は ★**非有界**（§2 `prod_of_prod_unbounded`）。

⇒ ★★**`Gained*` の値段を額面どおり 1 層分として積むと `AxLemma` は出ない。**
★`axDecay p k` **そのもの**（積ではなく単項）が 1 層の値段でなければならない。
★これは前波の `prod_depth_one_only_unbounded` と**同じ形の 2 例目**である。

## ★★正の結果 —— `AxWildDescent` の第 3 条件は**無料**（§3）

`AxTowerDecay.lean:473` の `AxWildDescent` は `x'` に 3 つを課す:

1. `wildDepth K x' < wildDepth K x`
2. `‖x − x'‖ ≤ c k · ε`
3. `∀ σ, ‖σ • x' − x'‖ ≤ c k · ε`   ←★**これ**

★★**3 は 1・2 と超距離＋等長だけから出る。** 実際

  `σ y − y = σ(y − x) + (σ x − x) + (x − y)`

で、`‖σ(y−x)‖ = ‖y−x‖`（等長 `AbsClosureModules.norm_smul_closure`）だから
超距離で 3 項の max、すなわち `max(‖x−y‖, ‖σx−x‖) ≤ c k · ε`。

⇒ ★★★**これから `AxWildDescent` を供給する者は「深さが下がる」と「近い」の 2 つだけ
示せばよい**（§4 `axWildDescent_of_stepwise_approx`）。
★木の `WildDepthDescent.exists_smul_invariant_of_padicValNat_index_succ`（`:517`）は
3 番目を**わざわざ証明している**（`norm_map_sum_quotient_sub_le` を使う 6 行）。
★★本補題があればその 6 行は要らない。★ただし既存ファイルは触らない規約なので置き換えていない。

## 逸脱の記録

- §3 の抽象核は `f : M → M` が加法的（差を保つ）かつ等長であることしか使わない。
  ★★**全単射性も乗法性も Galois も要らない。**
- ★本ファイルは `Gained*` の値段が「本当に `∏` なのか、それとも単項に落とせるのか」を
  **測っていない**。★§2 は「額面どおりに積むと閉じない」ことだけを言う。
- 前波のファイル `DepthApproxToAxLemma.lean` は**訂正しない**（他が読んでいる）。
  ★あちらの `axLemma_of_depth_bounds` は仮説付きの定理として**正しいまま**である。
  誤っていたのは「その仮説を `Gained*` が供給する」という**私の docstring の断定**であり、
  ★それを本ファイルで名指しして訂正した。
-/

namespace ABC3.Found.PGC

namespace StepwiseApproxFree

open ABC3.Skeleton.PGC

/-! ## §1 ★抽象核（その 1）—— 近似したあとの変位は増えない -/

section Kernel

variable {M : Type*} [SeminormedAddCommGroup M] [IsUltrametricDist M]

/-- ★★**近似したあとの変位は増えない**（超距離＋等長だけ）。

`f` が差を保ち等長で、`‖x − y‖ ≤ A` かつ `‖f x − x‖ ≤ A` なら `‖f y − y‖ ≤ A`。

★★`f` の**全単射性も乗法性も Galois も要らない**。
★これが `AxWildDescent` の第 3 条件（`ε` の伸び）が無料である理由である。 -/
theorem norm_sub_le_of_approx (f : M → M) (hsub : ∀ a b : M, f (a - b) = f a - f b)
    (hiso : ∀ a : M, ‖f a‖ = ‖a‖) {x y : M} {A : ℝ}
    (hxy : ‖x - y‖ ≤ A) (hx : ‖f x - x‖ ≤ A) : ‖f y - y‖ ≤ A := by
  have hrw : f y - y = f (y - x) + (f x - x) + (x - y) := by
    rw [hsub]; abel
  rw [hrw]
  refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ hxy)
  refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ hx)
  rw [hiso, norm_sub_rev]
  exact hxy

end Kernel

/-! ## §2 ★負の結果 —— 1 層に `∏` を払うと閉じない -/

section NotClosing

variable {p : ℕ} [Fact p.Prime]

/-- ★★★`Gained*` の値段 `∏_{j ∈ Icc 1 k} axDecay p j` を深さごとの 1 層分として積むと、
有限積は ★**非有界**。

⇒ `axLemma_of_wildDescent_Icc` の仮説（`∀ n, ∏_{Icc 1 n} c k ≤ C`）が満たせない。
★★1 層の値段は `axDecay p k` **そのもの**（積ではなく単項）でなければならない。

★前波の `WildStepFieldSupply.prod_depth_one_only_unbounded` と同じ形の 2 例目。 -/
theorem prod_of_prod_unbounded (C : ℝ) :
    ∃ n : ℕ, C < ∏ k ∈ Finset.Icc 1 n, (∏ j ∈ Finset.Icc 1 k, axDecay p j) := by
  have ha : (1 : ℝ) < axDecay p 1 := WildStepFieldSupply.one_lt_axDecay_one
  obtain ⟨n, hn⟩ := pow_unbounded_of_one_lt C ha
  refine ⟨n, lt_of_lt_of_le hn ?_⟩
  have hcard : (Finset.Icc 1 n).card = n := by simp
  calc axDecay p 1 ^ n = ∏ _k ∈ Finset.Icc 1 n, axDecay p 1 := by
        rw [Finset.prod_const, hcard]
    _ ≤ ∏ k ∈ Finset.Icc 1 n, ∏ j ∈ Finset.Icc 1 k, axDecay p j := by
        refine Finset.prod_le_prod (fun _ _ => le_of_lt (lt_trans zero_lt_one ha))
          (fun k hk => ?_)
        have hmem : (1 : ℕ) ∈ Finset.Icc 1 k := by
          simp only [Finset.mem_Icc] at hk ⊢
          omega
        have hsplit := Finset.mul_prod_erase (Finset.Icc 1 k) (fun j => axDecay p j) hmem
        have h1 : (1 : ℝ) ≤ ∏ j ∈ (Finset.Icc 1 k).erase 1, axDecay p j :=
          Finset.one_le_prod (fun i _ => one_le_axDecay p i)
        rw [← hsplit]
        exact le_mul_of_one_le_right (le_trans zero_le_one (one_le_axDecay p 1)) h1

end NotClosing

/-! ## §3 ★★正の結果 —— `AxWildDescent` の第 3 条件は無料 -/

section Free

variable {p : ℕ} [Fact p.Prime]

/-- §1 を `K.closure` の Galois 作用に代入したもの（等長は
`AbsClosureModules.norm_smul_closure`、`:270`）。 -/
theorem norm_smul_sub_le_of_approx (K : PAdicLocalField p) (σ : K.absGal)
    {x y : K.closure} {A : ℝ} (hxy : ‖x - y‖ ≤ A) (hx : ‖σ • x - x‖ ≤ A) :
    ‖σ • y - y‖ ≤ A :=
  norm_sub_le_of_approx (fun z => σ • z) (fun a b => smul_sub σ a b)
    (norm_smul_closure K σ) hxy hx

/-- ★★★**`AxWildDescent` を供給するには 2 つ示せばよい** —— 深さが下がることと、近いこと。

第 3 条件（`∀ σ, ‖σ • x' − x'‖ ≤ c k · ε`）は §1 から**自動**である。

★木の `WildDepthDescent.exists_smul_invariant_of_padicValNat_index_succ`（`:517`）は
第 3 条件を 6 行かけて証明している。★本補題があればそれは要らない
（規約により既存ファイルは触っていない）。 -/
theorem axWildDescent_of_stepwise_approx (K : PAdicLocalField p) {c : ℕ → ℝ}
    (hc : ∀ k, 1 ≤ c k)
    (hstep : ∀ (k : ℕ) (ε : ℝ), 0 ≤ ε → ∀ x : K.closure, wildDepth K x = k →
      (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
      ∃ x' : K.closure, wildDepth K x' < k ∧ ‖x - x'‖ ≤ c k * ε) :
    AxWildDescent K c := by
  intro ε hε x hx _hdvd
  obtain ⟨x', hlt, h1⟩ := hstep (wildDepth K x) ε hε x rfl hx
  exact ⟨x', hlt, h1, fun σ =>
    norm_smul_sub_le_of_approx K σ h1 ((hx σ).trans (le_mul_of_one_le_left hε (hc _)))⟩

end Free

/-! ## §4 使っている公理の一覧 -/

#print axioms norm_sub_le_of_approx
#print axioms prod_of_prod_unbounded
#print axioms norm_smul_sub_le_of_approx
#print axioms axWildDescent_of_stepwise_approx

end StepwiseApproxFree

end ABC3.Found.PGC
