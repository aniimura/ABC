import ABC3.Found.PGC.ZetaUnitFactor

/-!
# [pGC] ★★★★(Q2) は **#69 に触れずに落ちた** —— 柱の話が完全に閉じた

## ★★★まず自分の見積もりの訂正（本体の観察が正しかった）

本体が事実として提示したこと:

> ★数学の費用の見積もりは 6 波連続で当たっている。
> ★★**開いていないファイル・未着手の配管についての見積もりは 2 度続けて外した。**

★**3 度目を出さないため、本波は見積もらずに測った。**結果:

> (Q2)「`p ∤ b` なら `‖(b:F)‖ = 1`」は ★**Bézout と強三角不等式だけで落ちる**。
> `PAdicLocalField` も `adjoinIntegers` も**一切使わない**。15 行、**1 往復**、初回で通った。

⇒ ★★★**前波で「#69 の境界に近い」と書いたのは誤りだった。**
★配管についての見積もりを外したのは**これで 3 度目**である
（①「表記の変換 1 行」→ 中身のある主張、②「`v_L` への翻訳が要る」→ 要らなかった、
③「#69 の境界」→ 触れずに済んだ）。
★★**今後、開いていないものは「未測定」と書き、費用は書かない。**

## ★★核（分岐・付値・Galois の語が 1 語も出ない）

`norm_natCast_eq_one_ultra` : `F` が超距離ノルム体で `‖(p:F)‖ < 1`、`p ∤ b` なら
★**`‖(b:F)‖ = 1`**。

証明: `gcd(p,b) = 1` から Bézout `u·p + v·b = 1`（`Nat.isCoprime_iff_coprime`）。
`F` に送ると `‖1‖ ≤ max(‖u‖‖p‖, ‖v‖‖b‖)`。整数の像は `‖·‖ ≤ 1`
（`IsUltrametricDist.norm_intCast_le_one`）なので、もし `‖b‖ < 1` なら
右辺は `< 1` になって矛盾。★`‖b‖ ≤ 1` と合わせて `= 1`。

## ★★在庫の測定（#348 が効いた）

★最初 `norm_natCast_eq_one_of_not_dvd` という名前で書こうとしたら
`grep -rn "^theorem norm_natCast_eq_one_of_not_dvd" lean/ABC3/Found/PGC/*.lean` が **1 件**返した:

```
PrimeToPTorsion.lean:52  norm_natCast_eq_one_of_not_dvd (K : PAdicLocalField p) {m : ℕ}
                         (hm : ¬ (p ∣ m)) : ‖((m : ℕ) : K.carrier)‖ = 1
```

★★**木に既に在った** —— ただし `PAdicLocalField` 専用で、
`Padic.norm_natCast_eq_one_iff`（mathlib、`ℚ_p` 固有）を経由している。
★本波のものは**一般の超距離ノルム体**版なので**別物**だが、
短名が衝突するので `norm_natCast_eq_one_ultra` に改名した（#348）。
⇒ ★**「在る」と「使える」は別**の例がもう 1 つ増えた。

## ★★★これで柱の話は完全に閉じた

| 仮定 | 状態 |
|---|---|
| `ζ_{p^m} − 1` が素元 | ★定理（`ZetaSubOnePrime.norm_zeta_sub_one_pow`） |
| `v_L(μ) = p` | ★定理（`ZetaStepRatio.norm_zeta_step`） |
| `v_L(w) = p`、`v_L(σμ−μ) = p²` | ★定理（`ZetaUnitFactor`） |
| `‖(b:F)‖ = 1`（`p ∤ b`） | ★★**定理（本波）** |
| `‖(p:F)‖ < 1` | ★**設定の仮定**（`F` の剰余標数が `p`）。★穴ではない |

★★**残るのは `‖p‖ < 1` だけで、それは「どんな体で議論しているか」の宣言である。**
`[Fact p.Prime]` と同じ種類のもので、証明すべき命題ではない。
⇒ ★**本日の柱（`ρ = (1+π)w` の 2 成分性と `t = p(p−1)`）は仮定ゼロで閉じた。**

## ★残り（★見積もりは書かない。★測っていないものは「未測定」とだけ書く）

1. ①(a)「不分岐 ⇒ `∃ a, Tr(a) = 1`」 —— ★**測定済み、重い**
   （`UnramifiedStepFixed.lean:9-25`）。
2. ②′の一様な証人 / ③ / ⑤ —— ★変化なし。
3. `GainedTowerDescent.lean:94` の総当たりの表 —— ★**未測定**（開いていない）。
4. ★棚卸しは `ZetaUnitFactor.lean` の §0（前波）にある。本波は上の表 1 つを更新する。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
   ★`PrimeToPTorsion.lean:52` の同名宣言は**そのまま** —— 設定が違うので競合しない。
3. ★宣言名 2 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**で数え、★**1 件衝突していたので改名した**（#348 が実際に効いた）。
4. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace ResidueUnitNorm

/-! ## §1 抽象核 —— `p ∤ b` なら `‖b‖ = 1`（Bézout と強三角不等式だけ） -/

section Kernel

/-- ★★★**核 —— `p ∤ b` なら `‖(b:F)‖ = 1`。**

超距離ノルム体で `‖(p:F)‖ < 1` を仮定する。証明は **Bézout と強三角不等式だけ**で、
★`PAdicLocalField` も `adjoinIntegers` も使わない（★`lean-idioms.md` #69 に触れない）。

★木の `PrimeToPTorsion.lean:52 norm_natCast_eq_one_of_not_dvd` は
`PAdicLocalField` 専用（`Padic.norm_natCast_eq_one_iff` 経由）なので**別物**。
★短名が衝突するため改名した（#348）。 -/
theorem norm_natCast_eq_one_ultra {F : Type*} [NormedField F] [IsUltrametricDist F]
    {p b : ℕ} (hp : p.Prime) (hlt : ‖(p : F)‖ < 1) (hb : ¬ p ∣ b) :
    ‖(b : F)‖ = 1 := by
  obtain ⟨u, v, huv⟩ : IsCoprime (p : ℤ) (b : ℤ) :=
    Nat.isCoprime_iff_coprime.mpr ((Nat.Prime.coprime_iff_not_dvd hp).mpr hb)
  have hcast : (u : F) * (p : F) + (v : F) * (b : F) = 1 := by
    have h := congrArg (fun z : ℤ => (z : F)) huv
    push_cast at h
    exact h
  have hu : ‖(u : F)‖ ≤ 1 := IsUltrametricDist.norm_intCast_le_one F u
  have hv : ‖(v : F)‖ ≤ 1 := IsUltrametricDist.norm_intCast_le_one F v
  have hble : ‖(b : F)‖ ≤ 1 := IsUltrametricDist.norm_natCast_le_one F b
  by_contra hne
  have hblt : ‖(b : F)‖ < 1 := lt_of_le_of_ne hble hne
  have hmax : ‖(1 : F)‖ ≤ max ‖(u : F) * (p : F)‖ ‖(v : F) * (b : F)‖ := by
    rw [← hcast]
    exact IsUltrametricDist.norm_add_le_max _ _
  rw [norm_one] at hmax
  have h1 : ‖(u : F) * (p : F)‖ < 1 := by
    rw [norm_mul]
    nlinarith [norm_nonneg ((p : F)), norm_nonneg ((u : F))]
  have h2 : ‖(v : F) * (b : F)‖ < 1 := by
    rw [norm_mul]
    nlinarith [norm_nonneg ((b : F)), norm_nonneg ((v : F))]
  have := max_lt h1 h2
  linarith

end Kernel

/-! ## §2 帰結 —— `‖ξ^b − 1‖ = ‖ξ − 1‖` の仮定が `p ∤ b` だけになる -/

section Consequence

/-- ★★`‖ξ^b − 1‖ = ‖ξ − 1‖` の仮定が **`p ∤ b` だけ**になった。
⇒ 本日の柱（`ρ = (1+π)·w` の 2 成分性）が**仮定ゼロ**で閉じる。 -/
theorem norm_pow_sub_one_eq_of_not_dvd {F : Type*} [NormedField F] [IsUltrametricDist F]
    {p b : ℕ} (hp : p.Prime) (hlt : ‖(p : F)‖ < 1) (hb : ¬ p ∣ b)
    {ξ : F} (hξ : ‖ξ‖ ≤ 1) (hone : ‖ξ - 1‖ < 1) :
    ‖ξ ^ b - 1‖ = ‖ξ - 1‖ :=
  ZetaUnitFactor.norm_pow_sub_one_eq hξ (norm_natCast_eq_one_ultra hp hlt hb) hone

end Consequence

/-! ## §3 使っている公理の一覧 -/

#print axioms norm_natCast_eq_one_ultra
#print axioms norm_pow_sub_one_eq_of_not_dvd

end ResidueUnitNorm

end ABC3.Found.PGC
