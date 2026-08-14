/-!
# 器具の較正 — skeleton 方式が何を検査し、何を検査しないか

このファイルは常設の証拠であり、削除しない。`PLAN.md` §1「事実1」の実演そのもの。
mathlib 不要(素の Lean 4 のみ)なので、mathlib のビルド状態に関わらず常に検証できる。

## 要点

1. `theorem foo : P := by sorry` に対して Lean は**型検査しかしない**。
2. 論文の定義を1つ取り違えると、skeleton がビルドを通るだけでなく、
   後で `sorry` が「埋まって」しまい、**sorry無し・公理依存ゼロ**で
   `0 = 1` まで証明できる。
3. 検出手段は「仮説を満たす具体例を1つ構成せよ」(ゲート G2)。
4. 同じ誤りを `axiom` で書くと、この検出手段が存在しなくなる(ゲート G4)。

★ 原典の対象は存在する。空虚なのは**我々の転写**である。Lean はその差を報告しない。
-/

namespace ABC3.Meta.Calibration

/-! ## 1. 誤って転写した定義 -/

/-- 論文の定義を写した「つもり」の構造体。
    原文の条件を2つ写したが、この2つは同時には成り立たない。 -/
structure Tower where
  f          : Nat → Nat
  strictMono : ∀ n, f n < f (n + 1)
  bounded    : ∀ n, f n < 10

/-- 論文の Proposition を写した「つもり」の skeleton。
    ★ ビルドは通る。Lean は「型が付くか」しか見ていない。 -/
theorem prop_paper_skeleton (T : Tower) : T.f 0 = 7 := by
  sorry

/-! ## 2. この skeleton が何も言っていないことの証明 -/

theorem Tower.le_self (T : Tower) : ∀ n, n ≤ T.f n := by
  intro n
  induction n with
  | zero => exact Nat.zero_le _
  | succ k ih => exact Nat.succ_le_of_lt (Nat.lt_of_le_of_lt ih (T.strictMono k))

/-- `Tower` には実例が一つも存在しない。 -/
theorem no_tower : ¬ Nonempty Tower := by
  rintro ⟨T⟩
  exact absurd (T.bounded 10) (Nat.not_lt.mpr (T.le_self 10))

/-- ★ 同じ主張が `sorry` 無しで完全に証明できてしまう。 -/
theorem prop_paper_proved (T : Tower) : T.f 0 = 7 :=
  absurd (⟨T⟩ : Nonempty Tower) no_tower

/-- ★ 同じ資格で `0 = 1` も証明できる——
    つまり `prop_paper_proved` は論文について何も述べていない。 -/
theorem prop_paper_nonsense (T : Tower) : (0 : Nat) = 1 :=
  absurd (⟨T⟩ : Nonempty Tower) no_tower

-- 受理ゲート(ビルド通過・sorry無し・公理チェック)を全部通ることの実演:
-- どちらも "does not depend on any axioms" と表示される。
#print axioms prop_paper_proved
#print axioms prop_paper_nonsense

/-! ## 3. 検出する検査 — ゲート G2(非空虚 witness)

「仮説を満たす具体例を1つ構成せよ」と要求すれば、上の欠陥は skeleton の
時点で落ちる。`Tower` については構成が不可能であることを `no_tower` が
機械的に証明している。以下は、ゲートが「通る」側の形。 -/

/-- 修正版: 誤って写した有界性条件を落とす。 -/
structure Tower' where
  f          : Nat → Nat
  strictMono : ∀ n, f n < f (n + 1)

/-- 非空虚 witness。これが書けることが skeleton の受理条件(G2)。 -/
example : Nonempty Tower' :=
  ⟨{ f := id, strictMono := fun n => Nat.lt_succ_self n }⟩

/-! ## 4. なぜ `axiom` を禁じるのか — ゲート G4

未構築の基礎を `axiom` で埋めると、上の検出手段が存在しなくなる。 -/

/-- 未構築の基礎を `axiom` で埋めた場合(悪い例)。 -/
axiom bounded_ax : ∀ n : Nat, n < 10

/-- 偽の公理からは何でも出る。しかも「空虚か」を問う相手がいない。 -/
theorem anything_from_axiom : (0 : Nat) = 1 :=
  absurd (bounded_ax 10) (by decide)

-- "depends on axioms: [bounded_ax]" と表示される——依存は分かるが、
-- `bounded_ax` が偽であることは教えてくれない。真偽の判定は人手に戻る。
#print axioms anything_from_axiom

/-!
対して `structure` 版では「`Nonempty` が示せるか」という、機械にかけられる
具体的な問いへ変わる。

> **`axiom` は検出不能な破滅、`structure` は検出可能な空虚。**

ゆえに本プロジェクトでは、未構築の基礎は必ず `ABC3/Interface/` の `structure`
として置く。`axiom` 宣言はこのファイル(較正デモ)を唯一の例外として禁止し、
`tools/check.mjs` が機械的に検査する。
-/

end ABC3.Meta.Calibration
