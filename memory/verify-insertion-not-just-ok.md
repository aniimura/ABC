---
name: verify-insertion-not-just-ok
description: ファイルに宣言を足したら grep -c で数を確認する。Python の str.replace は黙って失敗し、検査器具は変更前のファイルを ok と言う
metadata:
  type: feedback
---

`ABC3` の Lean ファイルに宣言を追記するとき、**Python の `str.replace` が黙って失敗する**
(2026-08-18 に 4 回発生)。マーカー文字列が一致しないと無変化のまま `print('ok')` が出る。

★その直後に `node tools/leanfile.mjs <file>` を走らせると、**変更前のファイル**を検査して
`ok` を返す。★★結果として「一発で通った」と誤認し、記録・コミット・push まで進んだ。

**Why:** 検証器具が「変更が入ったこと」を確認していないため、
「無変更 + 検査 ok」と「変更 + 検査 ok」が区別できない。

**How to apply:**
- 追記は `sed -i 'Nr <file>'`(行番号指定)か `cat > <file> <<'EOF'`(全書き換え)を使う。
  どちらも確実に効く。Python の `str.replace` は使わない。
- 追記後は必ず `grep -c "^theorem \|^def "` などで**宣言数を数えて**入ったことを確認してから
  `leanfile.mjs` / `lake build` を走らせる。
- 「ok」だけを根拠に「通った」と書かない。

関連: [[heredoc-eats-backslashes]]
