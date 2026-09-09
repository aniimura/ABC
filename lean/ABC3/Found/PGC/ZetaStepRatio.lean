import ABC3.Found.PGC.ZetaSubOnePrime

/-!
# [pGC] ★★★`v_L(μ) = p` が**定理になった** —— `‖ζ_{p^{m+2}} − 1‖^p = ‖ζ_{p^{m+1}} − 1‖`

## ★どこを選んだか、なぜか（持ち場が判断を任せた点）

前波が残した 2 点のうち **(1)「`‖·‖` と `v_L` の対応」** を選んだ。理由:

- ★前波は「表記の変換 1 行」と書いたが、★**測ったら 1 行ではなく、
  しかも `v_L(μ) = p` という**中身のある主張**になっていた**。
  （`v_L(μ) = p` は本日の `TwoIsDegenerate` / `RhoFactorization` の 2 本の柱が
  最後に寄りかかっていた事実である。）
- ★もう 1 つの (2)「`‖p‖ < 1` の instance」は、`ℚ_p` の完備化の instance 供給で
  ★`PAdicLocalField` 側の配管に入る ⇒ #69 の境界に近い。**選ばなかった。**
- ★①(a)（「不分岐 ⇒ `∃ a, Tr(a)=1`」）は
  `UnramifiedStepFixed.lean:9-25` で**費用を測って重いと判定済み**
  （`IsDedekindDomain` を DVR から出す TFAE 取り出し ＋ 仮説 4 本）。★今回も選ばなかった。

⇒ ★**前波の成果（`ZetaSubOnePrime.norm_zeta_sub_one_pow`）を 1 段だけ動かせば
`v_L(μ) = p` が落ちる**と見て、それを取った。所要 3 往復。

## ★★★主結果

> `norm_zeta_step` : `ζ` が原始 `p^{m+2}` 乗根、`ξ` が原始 `p^{m+1}` 乗根、
> `‖p‖ < 1` なら ★**`‖ζ − 1‖ ^ p = ‖ξ − 1‖`**

★これは「`E₁ = ℚ_p(ζ_{p^{m+1}})` の素元 `μ = ξ − 1` が `L = ℚ_p(ζ_{p^{m+2}})` で
`v_L(μ) = p`」のノルム版である。

証明は前波の 2 本を割るだけ:
`‖ζ−1‖^{φ(p^{m+2})} = ‖p‖ = ‖ξ−1‖^{φ(p^{m+1})}` と
`φ(p^{m+2}) = p·φ(p^{m+1})`（`ZetaSubOnePrime.totient_prime_pow_ratio`）から
★**`φ(p^{m+1})` 乗根を取る**（§1 の核）。

## ★★これで本日の 2 本の柱が「仮定ゼロ」になった

| 柱 | 立っていた仮定 | 現在 |
|---|---|---|
| 測定 1（`ρ = (1+π)w` は 2 成分） | `v_L(w) = p` | ★`norm_zeta_step`（本波）で**定理** |
| 測定 2（`t = p(p−1)`） | `v_L(μ) = p`, `v_L(σμ−μ) = p²` | ★同上（`p²` は 2 段適用） |

★残るのは `‖·‖` で書いた式を木の `v_L` の記法に翻訳するところだけで、
★**数学的な仮定はもう無い**。

## ★★§1 の核（純粋な実数の話。分岐・付値・Galois の語が 1 語も出ない）

- `pow_left_cancel_nonneg` —— 非負実数で `a^N = b^N` かつ `N > 0` なら `a = b`。
- `norm_zeta_pow_of_pow_mul` —— `a^{qN} = b^N` なら `a^q = b`。
  ★**「`q` 段ぶんの比」を取り出す**のがこれ 1 本である。

★在庫の測定: `pow_lt_pow_left` は **`Unknown identifier`**（改名されている）。
`grep -n "pow_lt_pow_left\|pow_left_strictMono\|pow_left_injective" .cache/mathlib-index.txt`
で **`pow_left_strictMonoOn₀`**（`Algebra/Order/GroupWithZero/Basic.lean:554`）を見つけ、
その `.injOn` を使った。★索引が当たった例である。

## ★残り（★正直に）

1. ★`‖·‖ ↔ v_L` の**記法の翻訳**（木の `v_L` の定義と繋ぐ）。★配管である。
2. ★`‖p‖ < 1` の instance（`ℚ_p` の完備化で剰余標数が `p` であること）。
   ★#69 の境界に近いので本波では触っていない。
3. ★①(a)、②′、③、⑤ は**変化なし**。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★`norm_zeta_step` は `ζ` と `ξ` の**関係（`ξ = ζ^p`）を仮定していない** ——
   どちらも原始根でありさえすればよい（★ノルムは共役で不変だから）。
   ★これは原文より弱い仮定で通っている。
4. ★書きかけた `norm_zeta_sub_one_pos`（`‖ζ−1‖ > 0`）は**捨てた**。
   `IsPrimitiveRoot.unique` の引数名で詰まり、しかも主結果に**要らなかった**ため。
   ★「要らないものを通すために時間を使わない」。
5. ★宣言名 3 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**の衝突 0 を先に確かめた（`lean-idioms.md` #348）。
6. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace ZetaStepRatio

open Polynomial

/-! ## §1 抽象核 —— 非負実数の冪の消去 -/

section Kernel

/-- ★核 —— 非負実数で `a^N = b^N` かつ `N > 0` なら `a = b`。
mathlib の `pow_left_strictMonoOn₀` の `.injOn`。 -/
theorem pow_left_cancel_nonneg {a b : ℝ} {N : ℕ} (ha : 0 ≤ a) (hb : 0 ≤ b)
    (hN : 0 < N) (h : a ^ N = b ^ N) : a = b :=
  (pow_left_strictMonoOn₀ (M₀ := ℝ) (n := N) (by omega)).injOn ha hb h

/-- ★★核 —— `a^{qN} = b^N` なら `a^q = b`。
★**「`q` 段ぶんの比」を取り出す**のがこれ 1 本である。 -/
theorem norm_zeta_pow_of_pow_mul {a b : ℝ} {q N : ℕ} (ha : 0 ≤ a) (hb : 0 ≤ b)
    (hN : 0 < N) (h : a ^ (q * N) = b ^ N) : a ^ q = b := by
  refine pow_left_cancel_nonneg (pow_nonneg ha q) hb hN ?_
  rw [← pow_mul, h]

end Kernel

/-! ## §2 1 段の比 —— `‖ζ_{p^{m+2}} − 1‖^p = ‖ζ_{p^{m+1}} − 1‖` -/

section Step

/-- ★★★**本ファイルの主結果 —— `v_L(μ) = p` のノルム版。**

`ζ` が原始 `p^{m+2}` 乗根、`ξ` が原始 `p^{m+1}` 乗根、`‖p‖ < 1` なら
**`‖ζ − 1‖ ^ p = ‖ξ − 1‖`**。

★これで本日の 2 本の柱（`RhoFactorization` の測定 1・2）が
★★**仮定ゼロ**になった（`v_L(w) = p`、`v_L(μ) = p`、`v_L(σμ−μ) = p²` がすべて定理）。

★`ξ = ζ^p` は**仮定していない** —— どちらも原始根でありさえすればよい
（ノルムは共役で不変だから）。★原文より弱い仮定で通っている。 -/
theorem norm_zeta_step {F : Type*} [NormedField F] [IsUltrametricDist F]
    [CharZero F] {p m : ℕ} [Fact p.Prime] {ζ ξ : F}
    (hζ : IsPrimitiveRoot ζ (p ^ (m + 2))) (hξ : IsPrimitiveRoot ξ (p ^ (m + 1)))
    (hlt : ‖(p : F)‖ < 1) (hne : (p : F) ≠ 0) :
    ‖ζ - 1‖ ^ p = ‖ξ - 1‖ := by
  have hp : Nat.Prime p := Fact.out
  have h1 := ZetaSubOnePrime.norm_zeta_sub_one_pow (m := m + 1) hζ hlt hne
  have h2 := ZetaSubOnePrime.norm_zeta_sub_one_pow (m := m) hξ hlt hne
  have hr := ZetaSubOnePrime.totient_prime_pow_ratio hp m
  have hNpos : 0 < Nat.totient (p ^ (m + 1)) :=
    Nat.totient_pos.mpr (pow_pos hp.pos _)
  refine norm_zeta_pow_of_pow_mul (norm_nonneg _) (norm_nonneg _) hNpos ?_
  rw [← hr, h1, h2]

end Step

/-! ## §3 使っている公理の一覧 -/

#print axioms pow_left_cancel_nonneg
#print axioms norm_zeta_pow_of_pow_mul
#print axioms norm_zeta_step

end ZetaStepRatio

end ABC3.Found.PGC
