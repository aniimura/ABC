import ABC3.Skeleton.GenEll.Section2

/-!
# ★[GenEll] `Theorem 2.1` の載せ替えが取りこぼしたもの(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.11。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★★★★このファイルは「どこまで取れたか」の記録である

2026-08-27(第 426 ブロック)に `theorem_2_1` を構成へ載せ替え、`sorry` を外した。
☆そのとき取れたのは **(i) ⟹ (ii) だけ**であった。

★★★★**2026-09-02(第 1446)——両向きを取った**。
`Interface/GenEll/Thm21Setup.lean` の `Thm21Data.cover` 欄
（段 B の幾何「`X(ℚ̄)^{≤d}` の中の列は部分列を取れば
ある compactly bounded subset に丸ごと入る」）を仮定すれば、
`(ii) ⟹ (i)` の組み立ては完全に証明できる（`Found/GenEll/Thm21Equiv.lean`、`sorry` 0）。

## ★★★★それでも「絶対的に証明された」ではない

`Thm21Data.cover` は **`Interface` の欄であり未構成である**。
☆本ファイルの `thm21_cover_still_open` はその事実を型で残す。
★これは `Theorem 3.8` / `Corollary 4.3` / `Corollary 4.4` が
`EllModuliData` の欄を `Interface` で受けているのとまったく同じ状況である。

## ★★結論が空虚でないこと

`theorem_2_1_nonvacuous` で、実際に (i) が成り立つ場合に (ii) が出ることを確かめる。
★`theorem_2_1_converse_nonvacuous` で**逆向きも空虚でない**ことを確かめる。
-/

namespace ABC3.Check.GenEll

open ABC3.Found.GenEll ABC3.Interface.GenEll

/-- ★**最も単純な `Thm21Data`** —— 全体が compactly bounded だと思う場合。

☆これは「欄が満たせる形をしている」ことの確認であって、
原文の compactly bounded subset の理論ではない。 -/
def trivialThm21Data : Thm21Data ℕ :=
  { degLe := Set.univ
    CB := fun _ => True
    cover := fun _ _ => ⟨Set.univ, id, strictMono_id, trivial, fun _ => trivial⟩ }

/-- ★**結論は空虚ではない** —— 全体で成り立つ場合に部分集合でも出る。 -/
theorem theorem_2_1_nonvacuous :
    ∀ KV : Set ℕ, trivialThm21Data.CB KV →
      BDle (fun x : ↑(KV ∩ trivialThm21Data.degLe) => (1 + (0:ℝ)) * ((0:ℝ) + 0))
           (fun x : ↑(KV ∩ trivialThm21Data.degLe) => (0 : ℝ)) :=
  (ABC3.Skeleton.GenEll.theorem_2_1 trivialThm21Data
    (fun _ : ℕ => (0 : ℝ)) (fun _ => 0) (fun _ => 0) 0).mp ⟨0, fun _ => by norm_num⟩

/-- ★★**逆向きも空虚ではない** —— (ii) から (i) が実際に出る。 -/
theorem theorem_2_1_converse_nonvacuous :
    BDle (fun x : ↑trivialThm21Data.degLe => (1 + (0:ℝ)) * ((0:ℝ) + 0))
         (fun x : ↑trivialThm21Data.degLe => (0 : ℝ)) :=
  (ABC3.Skeleton.GenEll.theorem_2_1 trivialThm21Data
    (fun _ : ℕ => (0 : ℝ)) (fun _ => 0) (fun _ => 0) 0).mpr
    (fun _ _ => ⟨0, fun _ => by norm_num⟩)

/-- ★★★★**まだ構成されていないのは `Thm21Data.cover` 欄 1 本である** ——
その事実を型で残す。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

★★★★★**これは主張ではない。**「何が残っているか」を機械が読める形にしただけである。
残っているのは:

| 要るもの | 状態 |
|---|---|
| noncritical Belyi maps(`[NCBelyi] Theorem 2.5`) | 一般の曲線への帰着(Riemann–Roch)が未了 |
| `Proposition 1.7` | 局所から大域への組み立てが未了 |
| 双曲曲線の étale 基本群の構造 | `[Stacks] 58.6` に原典、大きさ未測定 |
| compactly bounded subset の support の理論 | `Example 1.3, (ii)`、posit のまま |

☆コンパクト性の側（`Found/GenEll/Thm21Extract.lean` の `exists_bad_seq_tendsto`）は
既に `sorry` 0 で手元にある。 -/
def thm21_cover_still_open : String :=
  "Thm21Data.cover 欄が未構成。noncritical Belyi maps・Proposition 1.7・étale 基本群の構造が要る"

end ABC3.Check.GenEll
