import ABC3.Found.FrdI.Def24

/-!
# [FrdI] Proposition 4.1 —— Primary Steps

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.75–p.77。

原文 (FrdI p.75):
> Proposition 4.1. (Primary Steps) Suppose further that C is of perfect and

## ★★この命題は「圏の条件をモノイドの条件に翻訳する」ものである

原文 (FrdI p.76):
> Proof. First, we consider assertion (i). By applying the second equivalence of

原文 (FrdI p.76):
> categories of Definition 1.3, (iii), (d), to the various pre-steps over A, it follows

★★原文の証明は毎条とも同じ形をしている:
**`Definition 1.3, (iii), (d)`(co-angular pre-step のコスライスと `Order(Φ(A))` の圏同値)を
当てて、圏の条件をモノイドの条件に書き換える**。そのうえでモノイドの側を解く。

★したがって**実装も 2 層に分ける**:

| 層 | 内容 |
|---|---|
| ★モノイド層 | 翻訳後の条件とモノイドの語(`primary` など)の同値。**本ファイル** |
| 圏層 | `coaPreUnderEquiv` による翻訳。`Definition 1.3, (iii), (d)` を要する |

★**モノイド層は §1–§2 の語彙だけで閉じる**ので、先にそこを取る。

## ★(i) のモノイド層(原文が明示している)

原文 (FrdI p.76):
> Proof. First, we consider assertion (i). By applying the second equivalence of

原文 (FrdI p.76):
> categories of Definition 1.3, (iii), (d), to the various pre-steps over A, it follows

★原文が「follows immediately」と言う中身は、★★**`⟸` で `perfect` を使って
`n` で割る**ところである。原文は「together with the fact that `Φ(A)` is perfect」と
一言添えるだけだが、**そこが唯一の非自明な段**である。
-/

namespace ABC3.Found.FrdI

universe w

variable {M : Type w} [AddCommMonoid M]

/-- ★**`n ≥ 1` なら `a ⪯ n • a`** —— `a + (n-1) • a = n • a`。 -/
theorem mle_nsmul_self_one {a : M} {n : ℕ} (hn : 0 < n) : MLe a (n • a) := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  exact ⟨m • a, (succ_nsmul' a m).symm⟩

/-- ★★★**[FrdI] Proposition 4.1, (i) のモノイド層** ——
`perfect` なモノイドでは、`x ≠ 0` について

  `x` が **primary** ⟺ `x = a + b`(`a, b ≠ 0`)なるどの分解でも `x ⪯ a`

原文 (FrdI p.76):
> Proof. First, we consider assertion (i). By applying the second equivalence of

★★**`⟹` は `perfect` を使わない**: `x = a + b` から `a ⪯ x` はすぐ出るので、
primary の定義をそのまま当てるだけ。

★★★**`⟸` が本体**である。`b ⪯ x`(すなわち `b + c = n • x`)から出発し、
**`perfect` で `n` 割り**して `b = n • b'`、`c = n • c'` と書くと `b' + c' = x`。
- `c' = 0` なら `b = n • x` なので `x ⪯ b` は直接出る
- `c' ≠ 0` なら仮定が `x ⪯ b'` を与え、`b' ⪯ b`(`b = n • b'`)と繋いで `x ⪯ b`

★**原文が「together with the fact that `Φ(A)` is perfect」と一言添える所**が、
この「`n` 割り」である。 -/
theorem isPrimaryElt_iff_of_perfect (hperf : IsPerfectMonoid M) {x : M} (hx : x ≠ 0) :
    IsPrimaryElt x ↔ ∀ a b : M, x = a + b → a ≠ 0 → b ≠ 0 → MPrec x a := by
  constructor
  · rintro ⟨-, hp⟩ a b hab ha -
    exact hp a ha ⟨1, one_pos, b, by rw [one_smul, ← hab]⟩
  · intro hcond
    refine ⟨hx, fun b hb hprec => ?_⟩
    obtain ⟨n, hn, c, hc⟩ := hprec
    -- ★`perfect` で `n` 割りする
    obtain ⟨b', hb'0⟩ := (hperf ⟨n, hn⟩).2 b
    obtain ⟨c', hc'0⟩ := (hperf ⟨n, hn⟩).2 c
    have hb' : n • b' = b := hb'0
    have hc' : n • c' = c := hc'0
    have hsum : (fun a : M => ((⟨n, hn⟩ : ℕ+) : ℕ) • a) (b' + c')
        = (fun a : M => ((⟨n, hn⟩ : ℕ+) : ℕ) • a) x := by
      show n • (b' + c') = n • x
      rw [smul_add, hb', hc', hc]
    have hx' : b' + c' = x := (hperf ⟨n, hn⟩).1 hsum
    have hbb' : MLe b' b := by
      rw [← hb']
      exact mle_nsmul_self_one hn
    by_cases hc0 : c' = 0
    · -- ★`c' = 0` —— `b' = x` なので `x ⪯ b` は直接
      have hbx : b' = x := by rw [← hx', hc0, add_zero]
      exact ⟨1, one_pos, by rw [one_smul, ← hbx]; exact hbb'⟩
    · -- ★`c' ≠ 0` —— 仮定を当てて繋ぐ
      have hb0 : b' ≠ 0 := by
        intro h0
        exact hb (by rw [← hb', h0, smul_zero])
      exact mprec_trans (hcond b' c' hx'.symm hb0 hc0) ⟨1, one_pos, by rw [one_smul]; exact hbb'⟩

/-! ### ★§4 を構造化した(2026-08-17)

★★**着手時点では `1_Structured/` に §4 の構造化 HTML が無かった**——
`frdi-thm-3-4` までしか sectionId が存在せず、ゲートの G1(locator 検証)を
通す `.src` を**一切書けなかった**。★これは数学の穴ではなく**器具の準備**の問題で、
★★**§4–§6 が 0 件で止まっていた理由の一部がここにあった**。

★`section-4.html` を作り、§4 の必要 10 項目すべてに sectionId を用意した
(逐語は抽出テキストから ASCII 純粋な断片のみを splice、S1–S6 PASS)。
★これで §4 の実装に `.src` が付けられる。
-/

/-- ★**locator** —— `Proposition 4.1, (i)` の**モノイド層**。

★★**条つき**(命題全体の完全実装の主張ではない)。圏層
(`Definition 1.3, (iii), (d)` による翻訳)は未実装なので数には出ない。 -/
def isPrimaryElt_iff_of_perfect.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 76, item := "Proposition 4.1, (i) — モノイド層",
    sectionId := "frdi-prop-4-1" }

end ABC3.Found.FrdI
